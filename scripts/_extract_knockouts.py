#!/usr/bin/env python3

import hail as hl
import argparse
import pandas as pd
import random
import string
import sys
import os.path

from ukb_utils import hail_init
from ko_utils import io
from ko_utils import ko


class SplitArgs(argparse.Action):
    r"""Method for splitting input csv into a list"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))


def get_tid(length=5):
    r"""method for getting random ID string for alleles"""
    return ''.join(random.choices(string.ascii_uppercase + string.ascii_lowercase, 
                                  k=length))

def main(args):

    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    out_type = args.out_type
    chrom = int(args.chrom)
    only_vcf = args.only_vcf
    checkpoint = args.checkpoint
    aggr_method = args.aggr_method

    export_all_gts = args.export_all_gts
    csqs_category = args.csqs_category

    # import phased/unphased data
    hail_init.hail_bmrc_init('logs/hail/knockout.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')

    mt = io.import_table(input_path, input_type, calc_info = False)

    # subset to current csqs category
    mt = mt.filter_rows(hl.literal(set(csqs_category)).contains(mt.consequence_category))
    n_csqs = mt.count()[0]
    sys.stderr.write(f"Filtering to {n_csqs} variants that are {csqs_category}.")
   
    # aggregate knockouts by gene 
    gene_expr = mt.consequence.vep.worst_csq_by_gene_canonical.gene_id
    genes = ko.collect_phase_count_by_expr(mt, gene_expr)
    genes = ko.sum_gts_entries(genes)
    expr_pko = ko.calc_prob_ko(genes.hom_alt_n, genes.phased, genes.unphased)
    expr_ko = ko.annotate_knockout(genes.hom_alt_n, expr_pko, genes.phased)

    # calculate probability of being knocked out based on phased counts
    genes = genes.annotate_entries(pKO=expr_pko, knockout=expr_ko)

    # write out variants involved and vcf
    genes = genes.filter_entries(hl.is_defined(genes.knockout)).entries()
    genes = genes.transmute(
            gts = hl.delimit(genes.gts, ";"),
            varid = hl.delimit(genes.varid, ";")
            )
    genes.flatten().export(out_prefix + "_all.tsv.gz")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--chrom', default=None, help='Chromosome to be used') 
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--only_vcf', default=False, action='store_true', help='Only return VCF (less memory required when running)')
    parser.add_argument('--checkpoint', default=False, action='store_true', help='Checkpoint gene-aggregation matrix to avoid Spark Memory overflow errors') 
    parser.add_argument('--aggr_method', default="collect", help='How should the CH matrix be generated?')
    # filtering options
    parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='What categories should be subsetted to?')

    args = parser.parse_args()

    main(args)



