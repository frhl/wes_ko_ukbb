#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ko_utils import ko
from ukb_utils import hail_init
from ukb_utils import variants

def main(args):
    
    phased_path = args.phased_path
    phased_type = args.phased_type
    prephased_path = args.prephased_path
    prephased_type = args.prephased_type
    parents_path = args.parents_path
    parents_type = args.parents_type
    out_prefix = args.out_prefix
    out_type = args.out_type

    # set up
    hail_init.hail_bmrc_init_local('logs/hail/merge_chunks.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    parents = io.import_table(parents_path, parents_type, calc_info = False)
    phased = io.import_table(phased_path, phased_type, calc_info = False)
    # if PS not in entry field, we need to set the phased genotype to unphased
    # as there is no relevant information for phasing singletons.
    if prephased_path:
        prephased = io.import_table(prephased_path, prephased_type)
        phased = phased.annotate_entries(
                PS = prephased[phased.row_key, phased.col_key].PS)
        expr_ambigious_singleton = (ko.is_phased(phased.GT)) & (
                variants.get_mac_expr(phased) == 0) & (hl.is_defined(phased.PS))
        n_ambigious_sites = phased.aggregate_entries(hl.agg.sum(expr_ambigious_singleton))
        if n_ambigious_sites > 0:
            print(f"Unphasing {n_ambigious_sites} singleton genotypes with ambigious phasing.")
            phased = phased.transmute_entries(
                GT=hl.if_else(
                    expr_ambigious_singleton,
                    ko.unphase(phased.GT),
                    phased.GT,
                    )
                )

    
    # PQ/PS in all phased data except chromsome 21
    if "PQ" in list(phased.col):
        parents = parents.annotate_entries(
                   PQ = hl.missing('int32'),
                   PS = hl.missing('int32')
                )
       
    
    mt = parents.union_cols(phased)
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--phased_path', default=None, help='What is the input directory of the files?')
    parser.add_argument('--phased_type', default="vcf", help='What is the input directory of the files?')
    parser.add_argument('--parents_path', default=None, help='What extension does the file(s) end with?')
    parser.add_argument('--parents_type', default="vcf", help='What extension does the file(s) end with?')
    parser.add_argument('--prephased_path', default=None, help='Was pre-phasing using? in that case annotate with PS field.')
    parser.add_argument('--prephased_type', default="vcf", help='Prephasing file type input')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default="vcf", help='vcf/plink/mt')
    args = parser.parse_args()

    main(args)


