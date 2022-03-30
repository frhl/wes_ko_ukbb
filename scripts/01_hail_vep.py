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
    chrom = args.chrom
    input_path = args.input_path
    input_type = args.input_type
    vep_path = args.vep_path
    out_prefix = args.out_prefix

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip

    # get file and annotate
    mt = qc.get_table(input_path=input_path, input_type=input_type) # 12788
    mt = process_consequences(hl.vep(mt, "utils/configs/vep_final.json"))
    ht = mt.rows()

    # write out VEP hail table
    ht.write(out_prefix + "_vep.ht", overwrite=True)

    # write a table with variants and consequences which
    # we can access in other coding languages
    ht.



if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--vep_path', default=None, help='Path to dbNSFP annotations')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--chrom', default=None, help='Chromosome to be used')

    args = parser.parse_args()

    main(args)



