#!/usr/bin/env python3

import hail as hl
import argparse
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

def main(args):

    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    csqs_category = args.csqs_category
    subset_gene = args.subset_gene
    chunk_idx = args.chunk_idx
    genes_per_chunk = args.genes_per_chunk

    # import phased/unphased data
    hail_init.hail_bmrc_init('logs/hail/_extract_knockouts.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
    
    checkpoint_dir = out_prefix + "_checkpoint.mt"
    checkpoint_file = checkpoint_dir + "/_SUCCESS"

    if not os.path.isfile(checkpoint_file):
        mt = io.import_table(input_path, input_type, calc_info=False)
        # create list for subsetting
        if subset_gene:
            subset_gene = list(set(subset_gene))
            subset_gene = list(filter(None, subset_gene))
            assert len(set(subset_gene)) >= 1
            print(subset_gene)
            gene_id_expr = mt.consequence.vep.worst_csq_by_gene_canonical.gene_id
            csqs_expr = hl.literal(set(csqs_category)).contains(mt.consequence_category)
            gene_expr = hl.literal(subset_gene).contains(gene_id_expr)
            mt = mt.filter_rows((csqs_expr & gene_expr))
        # checkpoint and write
        mt = mt.repartition(32)
        mt = mt.checkpoint(out_prefix + "_checkpoint.mt", overwrite = True)
    else:
        mt = hl.read_matrix_table(checkpoint_dir)

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
            gts=hl.delimit(genes.gts, ";"),
            varid=hl.delimit(genes.varid, ";")
            )
    genes.flatten().export(out_prefix + ".tsv.gz")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='What categories should be subsetted to?')
    parser.add_argument('--subset_gene', default=None, action=SplitArgs, help='Subset gene?')
    parser.add_argument('--chunk_idx', default=None, help='what chunk should be subset to?')
    parser.add_argument('--genes_per_chunk', default=None, help='How many genes per chunk?')

    args = parser.parse_args()

    main(args)



