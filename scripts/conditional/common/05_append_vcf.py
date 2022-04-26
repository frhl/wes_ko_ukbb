#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import tables
from ko_utils import io

def main(args):
    
    ko_path = args.ko_path
    mk_path = args.mk_path
    ko_type = args.ko_type
    mk_type = args.mk_type
    out_prefix = args.out_prefix
    out_type = args.out_type
   
    hail_init.hail_bmrc_init('logs/hail/append_vcf.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
   
    mk = io.import_table(mk_path, mk_type, calc_info=False)
    ko = io.import_table(ko_path, ko_type, calc_info=False)
    ko = tables.order_cols(ko, mk)
    mt = io.rbind_matrix_tables(ko, mk)
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mk_path', default=None, help='')
    parser.add_argument('--mk_type', default=None, help='')
    parser.add_argument('--ko_path', default=None, help='')
    parser.add_argument('--ko_type', default=None, help='')
    parser.add_argument('--out_prefix', default=None, help='')
    parser.add_argument('--out_type', default=None, help='')
    args = parser.parse_args()

    main(args)

