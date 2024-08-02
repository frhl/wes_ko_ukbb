#!/usr/bin/env python3

import hail as hl
import argparse
import sys
import os

from ukb_utils import hail_init
from ko_utils import io
from ko_utils import ko

class SplitArgs(argparse.Action):
    r"""Method for splitting input csv into a list"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))


def main(args):

    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    csqs_category = args.csqs_category

    # import phased/unphased data
    hail_init.hail_bmrc_init(log='logs/hail/knockout.log', default_reference='GRCh38', min_block_size=128)
    
    hl._set_flags(no_whole_stage_codegen='1')
    mt = io.import_table(input_path, input_type, calc_info = False)
    mt = mt.filter_rows(hl.literal(set(csqs_category)).contains(mt.consequence_category))
    n_csqs = mt.count()[0]
    gene_expr = mt.consequence.vep.worst_csq_by_gene_canonical.gene_id
    genes = ko.aggr_phase_count_by_expr(mt, gene_expr)
    genes = genes = genes.filter_entries(genes.phased.n > 0)
    genes.entries().flatten().export(out_prefix + "_alleles.txt.gz")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='What categories should be subsetted to?')

    args = parser.parse_args()

    main(args)



