#!/usr/bin/env python3

import hail as hl
import argparse

from gnomad.utils.vep import process_consequences
from ukb_utils import hail_init
from ko_utils import qc

def main(args):
    
    # parser
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    out_type   = args.out_type
    
    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/phasing_qc.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
  
    # Get phased data and reference
    mt = qc.get_table(input_path, input_type)
    rf = qc.get_table(input_reference, 'vcf')
    #rf = 
    
    



    # get tables
    mt = qc.get_table(input_path=input_path, input_type=input_type) 



if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_reference', default=None, help='Input to reference')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
 
    args = parser.parse_args()

    main(args)

