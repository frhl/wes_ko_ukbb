#!/usr/bin/env python3

import hail as hl
import argparse
import os

from ukb_utils import hail_init
from ko_utils import qc

def main(args):
    
    # parser
    input_path = args.input_path
    out_prefix = args.out_prefix
    
    hail_init.hail_bmrc_init('logs/hail/mt_to_vcf.log', 'GRCh38')
    mt = qc.get_table(input_path, 'mt')
    qc.export_table(mt, out_prefix = out_prefix, out_type = 'vcf')

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='Input path')
    parser.add_argument('--out_prefix', default=None, help='Output prefix')
 
    args = parser.parse_args()

    main(args)

