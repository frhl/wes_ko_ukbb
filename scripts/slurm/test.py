#!/usr/bin/env python3

import hail as hl
import argparse
import pandas
import os

from gnomad.utils.vep import process_consequences
from ukb_utils import hail_init
from ko_utils import qc
from ko_utils import analysis

def main(args):

    # parser
    arg1 = args.arg1
    arg2 = args.arg2
    print(arg1)
    print(arg2)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--arg1', default=None, help='Path to input')
    parser.add_argument('--arg2', default=None, help='Path to input')

    args = parser.parse_args()

    main(args)



