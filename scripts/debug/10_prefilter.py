#!/usr/bin/env python3

import hail as hl
import argparse
import pandas
import os

from gnomad.utils.vep import process_consequences
from ukb_utils import hail_init
from ko_utils import qc

def main(args):
    
    input_path = args.input_path
    input_type = args.input_type
    min_maf = args.min_maf
    out_prefix = args.out_prefix
    
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') #
    mt = qc.get_table(input_path=input_path, input_type=input_type) # 12788
    if min_maf is not None:
        mt = qc.filter_min_maf(mt, 0.001)
    hl.export_vcf(mt, out_prefix + '.vcf.bfz')

if __name__=='__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--min_maf', default=None, help='minimum MAF threshold') 
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


