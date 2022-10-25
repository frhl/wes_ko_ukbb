#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init

from ko_utils.variants import filter_min_mac, filter_missing

def main(args):

    # parser
    input_path = args.input_path
    input_type = args.input_type
    wes_sites = args.wes_sites
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/01_wes_gen.logs', 'GRCh38')
    mt = io.import_table(input_path, input_type, calc_info = False) # assuming build GRCh38
    if wes_sites:
        mt = mt.filter_rows(mt.wes == 1)
    io.export_table(mt, out_prefix, out_type) 

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='What is the input path to the file?')
    parser.add_argument('--input_type', default=None, help='What input type?')
    parser.add_argument('--wes_sites', default=None, action='store_true', help='Exclude parents of trio relationships')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type out vcf/plink/mt')
    args = parser.parse_args()

    main(args)


