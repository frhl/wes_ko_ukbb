#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import io

def main(args):
    
    in_file = args.in_file
    in_type = args.in_type
    out_prefix = args.out_prefix
    out_type = args.out_type
    extract_samples = args.extract_samples

    hail_init.hail_bmrc_init_local('logs/hail/extract_samples.log', 'GRCh38')
    samples = extract_samples.strip().split(",")
    mt = io.import_table(in_file, in_type)
    mt = mt.filter_cols(hl.literal(samples).contains(mt.s))
    io.export_table(mt, out_prefix, out_type)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--in_file', default=None, required = True, help='Path to input')
    parser.add_argument('--in_type', default=None, required = True, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--extract_samples', default=None, required = True, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, required = True, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, required = True, help='Type of output dataset (options: mt, vcf, plink)')
 
    args = parser.parse_args()

    main(args)

