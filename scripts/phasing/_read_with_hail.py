#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import variants
from ukb_utils import hail_init

def main(args):

    # parser
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    out_type = args.out_type
    calc_info = args.calc_info
    min_maf = args.min_maf

    hail_init.hail_bmrc_init_local('logs/hail/prephase_merge.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip

    mt=io.import_table(input_path, input_type, calc_info = False)
    if calc_info:
        mt=io.recalc_info(mt)
    if min_maf:
        mt = mt.filter_rows(variants.get_maf_expr(mt) >= float(min_maf))
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='Path prefix for intput dataset')
    parser.add_argument('--input_type', default=None, help='Either "mt", "vcf" or "plink"')
    parser.add_argument('--calc_info', default=None, action='store_true', help='Either "mt", "vcf" or "plink"')
    parser.add_argument('--min_maf', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Either "mt", "vcf" or "plink"')
    args = parser.parse_args()

    main(args)


