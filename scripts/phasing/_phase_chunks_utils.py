#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Nik Baya (2021-10-15)
"""

import hail as hl
import argparse
from datetime import datetime

from ukb_utils import RESOURCES_DIR
from ukb_utils.samples import get_fam_path
from ko_utils.io import recalc_info
#from ukb_wes_qc.utils import recalc_info

PHASE_WD = '/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb'
DATA_DIR = f'{PHASE_WD}/data'
TMP_DIR = f'{PHASE_WD}/data/tmp'


def get_prefilter_path_prefix(chrom):
    return f'{TMP_DIR}/prefilter/tmp-ukb_imp_phasing_prefilter_chr{chrom}'

def get_unphased_non_singleton_variants_path_prefix(chrom):
    return f'{TMP_DIR}/unphased/non_singl{chrom}'

def get_phasing_intervals_path(chrom, min_interval_unit):
    return f'{DATA_DIR}/intervals/intervals_min{min_interval_unit}_chr{chrom}.tsv'

def get_fam(app_id=11867, imp_200k_only=False):
    # use files created by the command:
    #   /well/lindgren/UKBIOBANK/nbaya/imp_200k/phase_ukb_wes/utils/make_fam.py --make_fam
    # github: https://github.com/nikbaya/phase_ukb_wes/blob/a5a97de8604d4ea1c7e7cfeecd9a1f57c2af2357/utils/make_fam.py
    fam_path = get_fam_path(
        app_id=app_id,
        imp_200k_only=imp_200k_only,
        relateds_only=False
    )
    fam = hl.import_table(
        paths=fam_path,
        key='IID',
        types={f: 'str' for f in ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHEN']}
    )
    return fam


def read_input(input_path, input_type, cache=False):
    assert input_type in {'mt', 'vcf', 'plink'}

    if input_type == 'mt':
        mt = hl.read_matrix_table(input_path)
    elif input_type == 'vcf':
        mt = hl.import_vcf(input_path, force_bgz=True,
                           array_elements_required=False)
    elif input_type == 'plink':
        mt = hl.import_plink(
            *[f'{input_path}.{x}' for x in ['bed', 'bim', 'fam']])

    if input_type != 'mt' and cache:
        mt = mt.cache()

    return mt


def maf_mac_filter(mt, maf=None, mac=None):
    '''Filter to variants with minor allele frequency > `maf` and
    minor allele count > `mac`
    For instance, use `mac` = 1 to remove singletons
    '''
    if mac is not None:
        if mt.info.AC.dtype == hl.dtype('array<int32>'):
            mt = mt.filter_rows(hl.min(mt.info.AC) > mac)
        elif mt.info.AC.dtype == hl.dtype('int32'):
            mt = mt.filter_rows((mt.info.AC > mac) & (
                mt.info.AC < mt.info.AN-mac))
        else:
            raise ValueError(
                'MatrixTable does not have an info.AC field with dtype int32 or array<int32>')
    if maf is not None:
        if mt.info.AF.dtype == hl.dtype('array<float64>'):
            mt = mt.filter_rows(hl.min(mt.info.AF) > maf)
        elif mt.info.AF.dtype == hl.dtype('float64'):
            mt = mt.filter_rows((mt.info.AF > maf) & (mt.info.AF < 1-maf))
        else:
            raise ValueError(
                'MatrixTable does not have an info.AF field with dtype float64 or array<float64>')
    return mt


def print_count(mt, out_prefix):
    '''Print output path prefix, variant count, sample count, and time
    '''
    n_variants, n_samples = mt.count()
    time = datetime.now().astimezone().strftime('%a %b %d %H:%M:%S %Z %Y')
    print(f'{out_prefix}: {n_variants} variants, {n_samples} samples\t{time}')


def count_variants(chrom, mt, maf=None, mac=None, printout=None):
    '''Count variants  that have minor allele frequency or minor allele count
    greater than `maf` or `mac`, respectively
    '''
    mt = maf_mac_filter(mt=mt,
                        maf=maf,
                        mac=mac)
    if printout is None:
        printout = f'chr{chrom} variants with maf>{maf}, mac>{mac}'
    print(printout+f': {mt.count_rows()}\n' +
          'path: '+get_filtered_path_prefix(chrom=chrom)+'.mt')


def filter_to_related(mt, get_unrelated=False, maf=None):
    '''Filter to samples in duos/trios, as defined by fam file
    '''
    fam = get_fam()
    # need to filter to samples present in mt first before collecting FIDs
    fam = fam.filter(hl.is_defined(mt.cols()[fam.IID]))
    fam_ct = fam.count()
    col_ct = mt.count_cols()
    # Samples with IIDs missing from the fam file cannot be included in the related sample set
    print(f'\n{col_ct-fam_ct}/{col_ct} samples have IIDs missing from fam file')
    iids = fam.IID.collect()
    # Get the FIDs of offspring with at least one parent in the dataset
    offspring_fids = fam.filter(
        hl.literal(iids).contains(fam.PAT)
        | hl.literal(iids).contains(fam.MAT)
    ).FID.collect()
    # Subset to samples which share FIDs with offspring
    fam = fam.filter(hl.literal(offspring_fids).contains(fam.FID))
    if get_unrelated:
        mt = mt.filter_cols(~hl.is_defined(fam[mt.s]))
    else:
        mt = mt.filter_cols(hl.is_defined(fam[mt.s]))
    mt = recalc_info(mt=mt, maf=maf)
    return mt


def annotate_with_fam_fields(mt):
    '''Annotate `mt` with FID (fam_id), PAT (pat_id), MAT (mat_id), SEX (is_female)
    from fam file with all UKB samples
    '''
    fam = get_fam(wes_200k_only=False)
    mt = mt.annotate_cols(fam_id=fam[mt.s].FID,
                          pat_id=fam[mt.s].PAT,
                          mat_id=fam[mt.s].MAT,
                          is_female=hl.if_else(fam[mt.s].SEX == '2',
                                               True,
                                               hl.if_else(fam[mt.s].SEX == '1',
                                                          False,
                                                          hl.null(hl.tbool))))
    return mt


def write_sites(mt, out_prefix, keep_fields=[]):
    '''Write sites-only tsv
    `keep_fields` is a list of fields to keep in addition to chrom, pos, ref, alt
    '''
    ht = mt.rows().key_by()
    ht = ht.select(chrom=ht.locus.contig,
                   pos=ht.locus.position,
                   ref=ht.alleles[0],
                   alt=ht.alleles[1],
                   *keep_fields)
    ht.export(f'{out_prefix}.tsv')


def write_phasing_intervals():
    """Write table of phasing intervals
    """
    print(f'Writing intervals to {out_prefix}')
    mt = mt.add_row_index()
    mt = mt.filter_rows(
        # take every `interval_size`-th variant
        ((mt.row_idx % interval_size) == 0)
        # include the last entry (subtract 1 because row_index is zero-based index)
        | (mt.row_idx == (mt.count_rows()-1))
    )
    write_sites(mt, out_prefix, keep_fields=['row_idx'])


def main(args):
    chrom = args.chrom
    if chrom == 'X' and args.chrX_filter_version != '':
        chrom = f'{chrom}-{args.chrX_filter_version}'

    if args.print_prefilter_path_prefix:
        print(get_prefilter_path_prefix(chrom))
    elif args.print_phased_non_singleton_variants_path_prefix:
        print(get_phased_non_singleton_variants_path_prefix(chrom))
    elif args.print_phasing_intervals_path:
        print(get_phasing_intervals_path(
            chrom=args.chrom,
            min_interval_unit=args.min_interval_unit
        ))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='Chromosome to filter')
    parser.add_argument('--print_prefilter_path_prefix',
                        action='store_true', help='Prints path prefix for phasing prefilter output')
    parser.add_argument('--print_pre_mendel_error_path_prefix',
                        action='store_true', help='Prints PLINK bfile path prefix for files used as input to Mendelian error filter')
    parser.add_argument('--print_mendel_error_path_prefix',
                        action='store_true', help='Prints PLINK bfile path prefix for PLINK files filtered for Mendelian errors')
    parser.add_argument('--print_mendel_error_stats_path_prefix',
                        action='store_true', help='Prints path prefix for PLINK to output stats on Mendelian errors')
    parser.add_argument('--print_phased_non_singleton_variants_path_prefix',
                        action='store_true', help='Prints path to phased VCF of all (non-singleton) variants for all samples')
    parser.add_argument('--chrX_filter_version', default='',
                        help="Version of filtering to run on chrX")
    parser.add_argument('--min_interval_unit', default=None,
                        help='Number of variants within each interval')
    parser.add_argument('--print_phasing_intervals_path', action='store_true',
                        help="Prints path to phasing intervals file. Requires defined args for --chrom, --phasing_region_size, --phasing_region_overlap and --max_phasing_region_size.")

    args = parser.parse_args()

    main(args)
