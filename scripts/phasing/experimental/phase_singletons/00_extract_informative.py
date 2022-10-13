#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import samples
from ukb_utils import variants
from ko_utils import io

def main(args):
    
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/annotate_mt.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')

    
    mt = io.import_table(in_file, in_type)
    mt = mt.annotate_rows(mac=variants.get_mac_expr(mt)) 
    print(mt.describe())



    # always produce a matrix-table 
    #io.export_table(mt, out_prefix, out_type)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--input_path', default=None, required = True, help='Path to input')
    parser.add_argument('--input_type', default=None, required = True, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, required = False, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, required = False, help='Type of output dataset (options: mt, vcf, plink)')
 
    args = parser.parse_args()

    main(args)

