#!/usr/bin/env python3

import hail as hl
import argparse
import pandas
import os

from gnomad.utils.vep import process_consequences
from ukb_utils import hail_init
from ko_utils import io

def main(args):

    # parser
    input_path = args.input_path
    out_prefix = args.out_prefix

    # standard revel usage
    ht = hl.read_table(input_path)
    ht = process_consequences(ht)
    ht.write(out_prefix + ".ht", overwrite=True)
    ht.export(out_prefix + ".txt.gz")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--out_prefix', default=None,help='Path prefix for output dataset')

    args = parser.parse_args()

    main(args)

