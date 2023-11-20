#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import samples


def main(args):

    # parser
    in_file = args.in_file
    in_type = args.in_type
    out_prefix = args.out_prefix

    hail_init.hail_bmrc_init_local('logs/hail/get_unrel.log', 'GRCh38')
    mt = io.import_table(in_file, in_type, calc_info=False)
    mt = samples.filter_ukb_to_unrelated_using_kinship(mt)
    mt.cols().export(out_prefix + ".txt.gz")


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file', default=None, help='Path for input dataset')
    parser.add_argument('--in_type', default=None, help='Type for file')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


