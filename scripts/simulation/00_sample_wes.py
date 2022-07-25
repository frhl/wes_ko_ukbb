#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import variants
from ukb_utils import samples


def main(args):

    # parser
    extract_samples = args.extract_samples
    exclude_samples = args.exclude_samples
    ancestry = args.ancestry
    exclude_related = args.exclude_related
    random_samples = args.random_samples
    random_seed = args.random_seed
    out_prefix = args.out_prefix
    out_type = args.out_type
    filter_to_unrelated_using_kinship_coef = args.filter_to_unrelated_using_kinship_coef
    filter_missing = args.filter_missing
    in_prefix = args.in_prefix
    in_type = args.in_type
    
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    mt = io.import_table(in_prefix, in_type)

    if extract_samples:
        ht_samples = hl.import_table(extract_samples, no_header=True, key='f0', delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_samples[mt.col_key]))
    
    if exclude_samples:
        ht_samples = hl.import_table(exclude_samples, no_header=True, key='f0', delimiter=',')
        mt = mt.filter_cols(~hl.is_defined(ht_samples[mt.col_key]))
  
    if exclude_related:
        related = samples.get_ukb_is_related_using_kinship_expr(mt)
        mt = mt.filter_cols(~related) 

    if filter_to_unrelated_using_kinship_coef:
        mt = samples.filter_ukb_to_unrelated_using_kinship(mt)

    if ancestry:
        mt = samples.filter_ukb_to_ancestry(mt, ancestry)

    if random_samples:
        mt = samples.choose_col_subset(mt, int(random_samples), seed = int(random_seed))

    if filter_missing:
        missing = hl.agg.mean(hl.is_missing(mt.GT)) <= float(filter_missing)
        mt = mt.filter_rows(missing)

    if out_type and out_prefix:
        io.export_table(mt, out_prefix, out_type)



if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--in_prefix', default=None, help='')
    parser.add_argument('--in_type', default=None, help='')
    parser.add_argument('--exclude_related', default=None, action='store_true', help='Exclude any related individuals.')
    parser.add_argument('--filter_to_unrelated_using_kinship_coef', default=None, action='store_true', help='Exclude any related individuals.')
    parser.add_argument('--filter_missing', default=None, help='Filter to variants with lt value in genotype missingness.')
    parser.add_argument('--extract_samples', default=None, help='Subset to sample IDs in MatrixTable')
    parser.add_argument('--exclude_samples', default=None, help='Exclude sample IDs from MatrixTable')
    parser.add_argument('--ancestry', default=None, help='Either "eur" or "all".')
    parser.add_argument('--random_samples', default=None, help='Subset to random samples')
    parser.add_argument('--random_seed', default=42, help='Seed for randomizer')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Either "mt", "vcf" or "plink"')
    args = parser.parse_args()

    main(args)


