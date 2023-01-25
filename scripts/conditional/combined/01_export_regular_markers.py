#!/usr/bin/env python3

import hail as hl
import argparse
import pandas as pd
import random
import string


from ukb_utils import hail_init
from ukb_utils import variants
from ko_utils import io
from ko_utils import ko


class SplitArgs(argparse.Action):
    r"""Method for splitting input csv into a list"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))


def main(args):

    mac_cutoff = args.mac_cutoff
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    out_type = args.out_type
    chrom = int(args.chrom)
    csqs_category = args.csqs_category

    # import phased/unphased data
    hail_init.hail_bmrc_init(log='logs/hail/knockout.log', default_reference='GRCh38', min_block_size=128)
    hl._set_flags(no_whole_stage_codegen='1')
    mt = io.import_table(input_path, input_type, calc_info = False)
    mt = mt.filter_rows(hl.literal(set(csqs_category)).contains(mt.consequence_category))
    mt = mt.annotate_rows(MAC = variants.get_mac_expr(mt))
    mt = mt.filter_rows(mt.MAC > int(mac_cutoff))
    mt = mt.annotate_entries(DS = mt.GT.n_alt_alleles())
    mt = mt.select_entries(mt.DS)
    if out_type not in "mt":
        mt = mt.checkpoint(out_prefix + ".mt", overwrite=True) 
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='Chromosome to be used') 
    parser.add_argument('--mac_cutoff', default=None, help='What mac_count should be used for collapsing') 
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='What categories should be subsetted to?')

    args = parser.parse_args()

    main(args)



