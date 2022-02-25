#!/usr/bin/env python3
"""
Prepares intervals for phasing non-singleton variants in chunks
Author: Nik Baya (2021-10-26)
"""

import argparse
import hail as hl
import pandas as pd

from ukb_utils import hail_init
from ukb_utils.vcf import import_vcf
#from phase_ukb_imputed.utils import get_unphased_non_singleton_variants_path_prefix, get_phasing_intervals_path

PHASE_WD = '/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb'
DATA_DIR = f'{PHASE_WD}/data'

#def get_phasing_intervals_path(chrom, min_interval_unit, interval_prefix):
#        return f'{interval_prefix}_min{min_interval_unit}.tsv'
#        #return f'{DATA_DIR}/intervals/intervals_min{min_interval_unit}_chr{chrom}.tsv'


#def get_unphased_non_singleton_variants_path_prefix(chrom):
#        return f'{DATA_DIR}/unphased/wes_union_calls/ukb_eur_wes_union_calls_200k_chr{chrom}'


def write_intervals(chrom: int, min_interval_unit: int, interval_path: str, vcf: str):
    """Write table with variants broken up into intervals
    Note that this is not used to get the intervals that we actually phase with,
    it's to get the chunks of variants that we then slice into to get intervals.
    :param chrom: Chromosome for which to write intervals
    :param min_interval_unit: (integer) Smallest units of variant counts in the interval
        files. E.g., if min_interval_unit=1000, the interval file would have every
        1000th variant printed out except for the final row, which has the last
        variant in the dataset.
    :param interval_prefix: path prefix to interval file
    """
    
    mt = import_vcf(
        path = vcf,
        reference_genome='GRCh38'
    )
    mt = mt.add_row_index()
    mt = mt.filter_rows(
        # Take every `min_interval_unit`-th variant and
        # include the last entry (subtract 1 because row_index is zero-based index)
        ((mt.row_idx % min_interval_unit) == 0)
        | (mt.row_idx == (mt.count_rows()-1))
    )
    ht = mt.rows().key_by()
    ht = ht.select(
        chrom=ht.locus.contig,
        pos=ht.locus.position,
        ref=ht.alleles[0],
        alt=ht.alleles[1],
        *['row_idx']
    )

    ht.export(interval_path)


def get_raw_start_idx(phasing_region_size: int, phasing_region_overlap: int,
                      phasing_idx: int):
    """Get the start_idx (left side of interval)
    This does not check the entire interval to see if it overlaps with the end
    of a chromosome.
    :param phasing_region_size: Number of variants in an interval to phase
    :param phasing_region_overlap: Overlap between phasing intervals in terms of
        variant count.
    :param phasing_idx: 1-based index for phasing intervals. The first interval
        has index = 1.
    """
    return int(phasing_region_size*(phasing_idx-1) - phasing_region_overlap*(phasing_idx-1))


def get_interval(chrom: str, min_interval_unit: int, phasing_region_size: int,
                 phasing_region_overlap: int, max_phasing_region_size: int, phasing_idx: int,
                 interval_path: str, reference_genome='GRCh38'):
    """Return genomic interval to be phased
    :param chrom: Chromosome of genomic interval
    :param min_interval_unit: (integer) Units of variant counts in the interval
        files. E.g., if min_interval_unit=1000, the interval file would have every
        1000th variant printed out.
    :param phasing_region_size: Number of variants in an interval to phase
    :param phasing_region_overlap: Overlap between phasing intervals in terms of
        variant count.
    :param max_phasing_region_size: Maximum number of variants to be included in
        a phasing interval.
    :param phasing_idx: 1-based index for phasing intervals. The first interval
        has index = 1.
    :param interval_prefix: path prefix to interval file
    :param reference_genome: Reference genome. Only affects whether a "chr"
        prefix is added to the region.
    """
    if reference_genome not in {'GRCh37', 'GRCh38'}:
        raise ValueError(f"Invalid reference genome: {reference_genome}")
    if max_phasing_region_size < phasing_region_size:
        raise ValueError("max_phasing_region_size must be greater than or equal to phasing_region_size")
    for arg_name, arg in [('phasing_region_size', phasing_region_size), ('phasing_region_overlap', phasing_region_overlap), ('max_phasing_region_size', max_phasing_region_size)]:
        if arg % min_interval_unit != 0:
            raise ValueError(f"{arg_name} must be a multiple of min_interval_unit")

    intervals = pd.read_csv(
        interval_path,
        #get_phasing_intervals_path(chrom, min_interval_unit, interval_prefix),
        sep='\t'
    )

    # 0-based index for variant count index (i.e. first variant will have index=0)
    start_idx = get_raw_start_idx(
        phasing_region_size=phasing_region_size,
        phasing_region_overlap=phasing_region_overlap,
        phasing_idx=phasing_idx
    )

    stop_idx = start_idx + phasing_region_size
    max_idx = intervals.row_idx.max()

    if start_idx > max_idx:
        return ''

    prev_phasing_idx = phasing_idx - 1
    prev_start_idx = get_raw_start_idx(
        phasing_region_size=phasing_region_size,
        phasing_region_overlap=phasing_region_overlap,
        phasing_idx=prev_phasing_idx
    )
    # If not the first interval and previous interval could be extended to max size
    if (phasing_idx > 1) and ((max_idx - prev_start_idx) <= max_phasing_region_size):
        return ''

    # Check if interval can be extended to max_phasing_region_size to cover the
    # final variants in a chromosome
    if (stop_idx < max_idx) and ((max_idx - start_idx) <= max_phasing_region_size):
        stop_idx = max_idx

    # If interval overhangs the end of the chromosome
    if stop_idx > max_idx:
        stop_idx = max_idx
        # Move the start_idx back so that we're still phasing an interval of size=phasing_region_size
        tmp_start_idx = stop_idx - phasing_region_size
        if tmp_start_idx < 0:
            start_idx = 0
        else:
            # Round down to nearest multiple of min_interval_unit
            start_idx -= start_idx % min_interval_unit

    start_pos = intervals.loc[intervals.row_idx == start_idx, 'pos'].values[0]
    stop_pos = intervals.loc[intervals.row_idx == stop_idx, 'pos'].values[0]

    if (phasing_region_overlap == 0) and (stop_idx != max_idx):
        stop_pos -= 1

    chrom_prefix = 'chr' if reference_genome == 'GRCh38' else ''

    return f"{chrom_prefix}{chrom}:{start_pos}-{stop_pos}"


def get_max_phasing_idx(chrom, min_interval_unit, phasing_region_size, phasing_region_overlap, max_phasing_region_size, interval_path):
    """Get maximum phasing index
    """
    # TODO: Make this more elegant + efficient
    phasing_idx = 1
    kwargs = {
        'chrom': chrom,
        'min_interval_unit': min_interval_unit,
        'phasing_region_size': phasing_region_size,
        'phasing_region_overlap': phasing_region_overlap,
        'max_phasing_region_size': max_phasing_region_size,
        'interval_path': interval_path
    }
    while get_interval(phasing_idx=phasing_idx, **kwargs) != '':
        phasing_idx += 1
        if phasing_idx > 100:
            break
    return phasing_idx - 1


def main(args):

    min_interval_unit = int(args.min_interval_unit)
    phasing_region_size = int(
        args.phasing_region_size) if args.phasing_region_size is not None else None
    phasing_region_overlap = int(
        args.phasing_region_overlap) if args.phasing_region_overlap is not None else None
    max_phasing_region_size = int(
        args.max_phasing_region_size) if args.max_phasing_region_size is not None else None
    phasing_idx = int(
        args.phasing_idx) if args.phasing_idx is not None else None

    if args.write_intervals:
        hail_init.hail_bmrc_init_local(
            log='logs/phase_chunks_hail.log',
            default_reference=args.reference_genome
        )

        write_intervals(
            chrom=args.chrom,
            min_interval_unit=min_interval_unit,
            interval_path=args.interval_path,
            vcf=args.target_vcf
        )

    elif args.get_interval:
        interval = get_interval(
            chrom=args.chrom,
            min_interval_unit=min_interval_unit,
            phasing_region_size=phasing_region_size,
            phasing_region_overlap=phasing_region_overlap,
            max_phasing_region_size=max_phasing_region_size,
            phasing_idx=phasing_idx,
            interval_path=args.interval_path,
            reference_genome=args.reference_genome
        )
        print(interval)
    elif args.get_max_phasing_idx:
        max_phasing_idx = get_max_phasing_idx(
            chrom=args.chrom,
            min_interval_unit=min_interval_unit,
            phasing_region_size=phasing_region_size,
            phasing_region_overlap=phasing_region_overlap,
            max_phasing_region_size=max_phasing_region_size,
            interval_path=args.interval_path
        )
        print(max_phasing_idx)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None,
                        help='Chromosome for which to write intervals')
    parser.add_argument('--min_interval_unit',
                        help='Number of variants within each interval. Interval file will have every min_interval_unit-th variant printed out.')
    parser.add_argument('--write_intervals', action='store_true',
                        help='Write file with intervals of variants, with each row a multiple of min_interval_unit')
    parser.add_argument('--phasing_region_size', default=None,
                        help="Desired number of variants to be included in an interval to phase")
    parser.add_argument('--phasing_region_overlap', default=None,
                        help="Overlap between phasing intervals in terms of variant count")
    parser.add_argument('--max_phasing_region_size', default=None,
                        help="Maximum number of variants to be included in a phasing interval.")
    parser.add_argument('--phasing_idx', default=None,
                        help="1-based index for phasing intervals. The first interval has index = 1")
    parser.add_argument('--reference_genome',
                        default='GRCh37', help="Reference genome")
    parser.add_argument('--get_interval', action='store_true',
                        help="Get genomic interval for phasing")
    parser.add_argument('--get_max_phasing_idx', action='store_true',
                        help="Get maximum phasing index")
    parser.add_argument('--interval_path', default=None,
                        help="Path prefix to interval file")
    parser.add_argument('--target_vcf', default=None,
                        help="Path prefix to vcf file for which to write intervals")



    args = parser.parse_args()

    main(args)

