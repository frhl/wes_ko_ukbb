#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init

def main(args):
    
    phased_path = args.phased_path
    phased_type = args.phased_type
    out_prefix = args.out_prefix
    hail_init.hail_bmrc_init_local('logs/hail/write_ps.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    mt = io.import_table(phased_path, phased_type, calc_info = False)
    mt.s.export(out_prefix + ".samples", header=False)

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--phased_path', default=None, help='What is the input directory of the files?')
    parser.add_argument('--phased_type', default="vcf", help='What is the input directory of the files?')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


