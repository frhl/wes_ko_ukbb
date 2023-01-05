#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import tables
from ko_utils import io

def main(args):
    
    chrom = args.chrom
    ko_path = args.ko_path
    ko_type = args.ko_type
    common_path = args.common_path
    common_type = args.common_type
    rare_path = args.rare_path
    rare_type = args.rare_type
    out_prefix = args.out_prefix
    out_type = args.out_type
   
    hail_init.hail_bmrc_init('logs/hail/combine_rare_common.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
 
    # there will always be rare variant
    ko = io.import_table(ko_path, ko_type, calc_info=False)
    rare = io.import_table(rare_path, rare_type, calc_info=False)
    rare = tables.order_cols(ko, rare)
    combined = io.rbind_matrix_tables(rare, ko)
    
    # there will not always be common variants
    if common_path:
        contig = "chr" + chrom
        common = io.import_table(common_path, common_type, calc_info=False)
        common = common.filter_rows(common.locus.contig == contig)
        common = tables.order_cols(ko, common)
        combined = io.rbind_matrix_tables(common, combined)

    # merge tables and export
    io.export_table(combined, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='')
    parser.add_argument('--ko_path', default=None, help='')
    parser.add_argument('--ko_type', default=None, help='')
    parser.add_argument('--common_path', default=None, help='')
    parser.add_argument('--common_type', default=None, help='')
    parser.add_argument('--rare_path', default=None, help='')
    parser.add_argument('--rare_type', default=None, help='')
    parser.add_argument('--out_prefix', default=None, help='')
    parser.add_argument('--out_type', default=None, help='')
    args = parser.parse_args()

    main(args)

