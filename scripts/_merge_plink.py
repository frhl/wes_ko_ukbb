#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init

def main(args):

    # parser
    in_prefix = args.in_prefix
    in_type = args.in_type
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip

    AUTOSOMES = list(map(str, range(1, 23)))
    mts = [io.import_table(in_prefix + "_chr" + str(chrom), in_type) for chrom in AUTOSOMES]
    mt = hl.MatrixTable.union_rows(*mts)
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--in_prefix', default=None, help='Path prefix for intput dataset')
    parser.add_argument('--in_type', default=None, help='Either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Either "mt", "vcf" or "plink"')
    args = parser.parse_args()

    main(args)


