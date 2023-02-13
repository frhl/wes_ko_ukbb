#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import samples

def main(args):
    
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    random_samples_count = args.random_samples_count
    seed = args.seed

    hail_init.hail_bmrc_init_local('logs/hail/random_samples.log', 'GRCh38')
    mt = io.import_table(input_path, input_type, calc_info = False)
    mt = samples.choose_col_subset(mt, int(random_samples_count), int(seed))    
    io.export_table(mt, out_prefix, "mt")


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='')
    parser.add_argument('--input_type', default=None, help='')
    parser.add_argument('--out_prefix', default=None, help='')
    parser.add_argument('--seed', default=None, help='')
    parser.add_argument('--random_samples_count', default=None, help='')

    args = parser.parse_args()

    main(args)


