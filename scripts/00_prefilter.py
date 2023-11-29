#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import samples
from ukb_utils import variants
from ko_utils import io

def main(args):
    # Initialize Hail and import data
    hail_init.hail_bmrc_init('logs/hail/knockout.log', 'GRCh38')
    mt = io.import_table(args.input_path, args.input_type, calc_info=False)

    print(mt.describe())

    # Apply filters based on the provided arguments
    if args.sex != 'both':
        mt = samples.filter_to_sex(mt, args.sex)

    # Filter based on final sample list
    if args.final_sample_list:
        ht_final_samples = hl.import_table(
            args.final_sample_list,
            no_header=True, key='f0',
            delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))

    # Filter based on final variant list
    if args.final_variant_list:
        ht_final_variants = hl.import_table(
            args.final_variant_list,
            types={'locus': hl.tlocus(reference_genome='GRCh38'),
                   'alleles': hl.tarray(hl.tstr)})
        ht_final_variants = ht_final_variants.key_by(
            ht_final_variants.locus, ht_final_variants.alleles)
        mt = mt.filter_rows(hl.is_defined(ht_final_variants[mt.row_key]))

    if args.maf_min:
        maf_min = float(args.maf_min)
        if maf_min > 0:
            maf_min_expr = variants.get_maf_expr(mt)
            mt = mt.filter_rows(maf_min_expr >= maf_min)

    if args.maf_max:
        maf_max = float(args.maf_max)
        if maf_max < 1:
            maf_max_expr = variants.get_maf_expr(mt)
            mt = mt.filter_rows(maf_max_expr <= maf_max)

    if args.exclude:
        ht = hl.import_table(args.exclude, impute=True).key_by('varid')
        mt = mt.filter_rows(~hl.literal(set(ht.varid.collect())).contains(mt.rsid))

    if args.pp_cutoff:
        pp_cutoff = float(args.pp_cutoff)
        expr_pp_cutoff = (mt.PP >= pp_cutoff) & (hl.is_defined(mt.GT))
        expr_keep = ~(hl.is_defined(mt.PP)) & (hl.is_defined(mt.GT))
        mt = mt.filter_entries(expr_pp_cutoff | expr_keep)

    # re-calculate info and export the processed data
    mt = io.recalc_info(mt)
    io.export_table(mt, args.out_prefix, args.out_type)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # Define command-line arguments
    parser.add_argument('--input_path', required=True, help='Path to input')
    parser.add_argument('--input_type', required=True, choices=['mt', 'vcf', 'plink'], help='Input type')
    parser.add_argument('--out_prefix', required=True, help='Path prefix for output dataset')
    parser.add_argument('--out_type', required=True, choices=['mt', 'vcf', 'plink'], help='Type of output dataset')
    parser.add_argument('--sex', default='both', choices=['male', 'female', 'both'], help='Filter to sex')
    parser.add_argument('--maf_min', type=float, help='Select all variants with a MAF greater than this value')
    parser.add_argument('--maf_max', type=float, help='Select all variants with a MAF less than this value')
    parser.add_argument('--exclude', help='Exclude variants by rsid and/or variant id')
    parser.add_argument('--pp_cutoff', type=float, help='PP cutoff value for filtering')
    parser.add_argument('--final_sample_list', help='Path to the final list of samples for filtering')
    parser.add_argument('--final_variant_list', help='Path to the final list of variants for filtering')


    args = parser.parse_args()
    main(args)

