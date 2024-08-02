#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import qc

def main(args):
    
    # parser
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    
    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    mt = qc.get_table(input_path=input_path, input_type=input_type)

    # remove invariant sites
    mt = hl.variant_qc(mt, name='variant_qc')
    mt = mt.filter_rows((mt.variant_qc.AF[0] > 0.0) & (mt.variant_qc.AF[0] < 1.0))

    # filter rows and write
    ht_rows_filter = mt.rows()
    ht_rows_filter.select().write(out_prefix, overwrite=True)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
 
    args = parser.parse_args()

    main(args)

