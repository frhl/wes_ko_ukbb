#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init

def main(args):

    out_prefix = args.out_prefix
    hail_init.hail_bmrc_init_local('logs/hail/annotate_mt.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--out_prefix', default=None, required = True, help='Path prefix for output dataset')
 
    args = parser.parse_args()

    main(args)

