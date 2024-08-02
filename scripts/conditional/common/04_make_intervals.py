#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init


def main(args):

    out_prefix = args.out_prefix
    gene_table = args.gene_table
    padding = int(float(args.padding))

    reference_genome = 'GRCh38'
    hail_init.hail_bmrc_init_local(
        'logs/hail/hail_conditional_tables.log',
        reference_genome)
    hl._set_flags(no_whole_stage_codegen='1')

    # import table of start/end contigs
    ht = hl.import_table(
        gene_table,
        force=True,
        delimiter='\t',
        missing='',
        types={
            'contig': hl.tstr,
            'start': hl.tint,
            'end': hl.tint})
    ht = ht.annotate(contig=hl.delimit(['chr', ht.contig], ''))
    
    # get refernece genome with lengths
    rg = hl.get_reference('GRCh38')
    chr_lens = hl.literal(rg.lengths)
    ht = ht.annotate(contig_length=chr_lens[ht.contig])

    # annotate regions with padding
    ht = ht.annotate(ranges=hl.array([ht.start, ht.end]))
    ht = ht.annotate(
        end=hl.max(ht.ranges),
        start=hl.min(ht.ranges))
    ht = ht.annotate(
        end_with_padding=hl.min(hl.array([ht.end + padding, ht.contig_length])),
        start_with_padding=hl.max(hl.array([ht.start - padding, 1])))
    ht = ht.drop(ht.ranges)

    # check if locus is valid, i.e. does it go above contig
    # length? This is only an issue when padding is large
    ht = ht.annotate(
        start_valid=hl.is_valid_locus(
            ht.contig,
            ht.start_with_padding,
            reference_genome))
    ht = ht.annotate(
        end_valid=hl.is_valid_locus(
            ht.contig,
            ht.end_with_padding,
            reference_genome))
    ht = ht.annotate(
            valid_intervals=ht.start_valid & ht.end_valid)

    defined_coords = hl.is_defined(ht.valid_intervals)
    n_drop = ht.filter(~defined_coords).count()
    if n_drop > 0:
        print(f'Dropping {n_drop} undefined start/end genomic coordinates.')
    ht = ht.filter(defined_coords)
    invalid_intervals = ht.filter(~ht.valid_intervals).count()
    assert invalid_intervals == 0, 'Some of the supplied intervals are outside current chromosome contig'
    
    # annotate intervals
    ht = ht.annotate(
        intervals=hl.locus_interval(
            ht.contig,
            ht.start_with_padding,
            ht.end_with_padding,
            True,
            True,
            reference_genome=reference_genome))

    # write hail talbe
    ht.write(out_prefix + ".ht")
    ht.export(out_prefix + ".txt.gz")

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--padding', default=0, help='How much extra padding should be included around genes (upstream/downsream)?')
    parser.add_argument('--gene_table', default=None, help='Path to HailTable that contains the significant genes from primary analysis.')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset (plink format)')
    args = parser.parse_args()

    main(args)
