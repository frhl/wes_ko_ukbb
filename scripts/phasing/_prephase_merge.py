#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init

def main(args):

    # parser
    input_list = args.input_list
    input_type = args.input_type
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/prephase_merge.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip

    # read list of files to be merged
    lines = []
    with open(input_list, "r") as infile:
        for line in infile:
            if line.strip():
                lines += [line.strip()]

    mts = [io.import_table(infile, input_type) for infile in lines]
    mt = hl.MatrixTable.union_cols(*mts)
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_list', default=None, help='Path prefix for intput dataset')
    parser.add_argument('--input_type', default=None, help='Either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Either "mt", "vcf" or "plink"')
    args = parser.parse_args()

    main(args)


