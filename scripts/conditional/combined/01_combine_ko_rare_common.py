#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import tables
from ko_utils import io

def main(args):
    
    common_path = args.common_path
    common_type = args.common_type
    ko_rare_path = args.ko_rare_path
    ko_rare_type = args.ko_rare_type
    out_prefix = args.out_prefix
    out_type = args.out_type
   
    hail_init.hail_bmrc_init('logs/hail/combine_ko_rare_common.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
  
    # assuming ko_rare contains psuedo-variants and actual variants
    ko_rare = io.import_table(ko_rare_path, ko_rare_type, calc_info=False)
    common = io.import_table(common_path, common_type, calc_info=False)
   
    # merge tables and export
    ko_rare = tables.order_cols(ko_rare, common)
    final = io.rbind_matrix_tables(ko_rare, common)
    io.export_table(final, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--common_path', default=None, help='')
    parser.add_argument('--common_type', default=None, help='')
    parser.add_argument('--ko_rare_path', default=None, help='')
    parser.add_argument('--ko_rare_type', default=None, help='')
    parser.add_argument('--out_prefix', default=None, help='')
    parser.add_argument('--out_type', default=None, help='')
    args = parser.parse_args()

    main(args)

