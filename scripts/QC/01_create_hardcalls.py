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
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
  
    # get tables
    mt = qc.get_table(input_path=input_path, input_type=input_type) 

    # Run VEP + gnomAD variant annotations
    #mt = process_consequences(hl.vep(mt, "utils/configs/vep_env.json"))
    
    # export files
    mt.select_entries(mt.GT).repartition(512).write(out_prefix + '.mt', overwrite=True)

    # counts
    n = mt.count()
    print('n samples:')
    print(n[1])
    print('n variants:')
    print(n[0])

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
 
    args = parser.parse_args()

    main(args)

