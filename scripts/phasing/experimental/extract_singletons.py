#!/usr/bin/env python3

import hail as hl
import argparse
import math
import os
import re

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import variants


def main(args):
    
    input_path = args.input_path
    input_type = args.input_type
    extract_samples = args.extract_samples
    out_prefix = args.out_prefix
    out_type = args.out_type

    # set up
    hail_init.hail_bmrc_init_local('logs/hail/extract_singletons.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
   
    # create list of samples
    extract_samples = extract_samples.strip().split(",")

    # filter to singeltons and low allele varaints
    mt = io.import_table(input_path, input_type, calc_info = False)
    mt = mt.filter_cols(hl.literal([extract_samples]).contains(mt.s))
    mt = mt.annotate_rows(mac=variants.get_mac_expr(mt))
    mt = mt.filter_rows(mt.mac<5)
    print(mt.describe())


    #io.export_table(mt, out, out_type)

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='What is the input directory of the files?')
    parser.add_argument('--input_type', default=None, help='What extension does the file(s) end with?')
    parser.add_argument('--extract_samples', default=None, help='What extension does the file(s) end with?')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='vcf/plink/mt')
    args = parser.parse_args()

    main(args)


