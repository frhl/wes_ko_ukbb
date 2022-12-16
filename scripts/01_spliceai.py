#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import io

def main(args):

    # parser
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/spliceai.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip

    # get file and annotate
    mt = io.import_table(input_path=input_path, input_type=input_type, calc_info = False) # 12788
    ht = mt.rows()
    ht = ht.drop(ht.info)
    hl.export_vcf(ht, out_prefix + ".vcf", tabix=True)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')

    args = parser.parse_args()

    main(args)

