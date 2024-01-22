#!/usr/bin/env python3

import hail as hl
import argparse
import pandas as pd
import random
import string
import sys
import os

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
    repartition = args.repartition
    export_all_gts = args.export_all_gts
    csqs_category = args.csqs_category
    exclude_singletons = args.exclude_singletons
    only_singletons = args.only_singletons
    print(args)

    # import phased/unphased data
    hail_init.hail_bmrc_init(log='logs/hail/knockout.log', default_reference='GRCh38', min_block_size=128)
    hl._set_flags(no_whole_stage_codegen='1')
    
    sys.stderr.write("No checkpoint found. Reading raw genotypes..")
    mt = io.import_table(input_path, input_type, calc_info = False)
    
    # subset to current csqs category
    mt = mt.filter_rows(hl.literal(set(csqs_category)).contains(mt.consequence_category))
    n_csqs = mt.count()[0]
    sys.stderr.write(f"Filtering to {n_csqs} variants that are {csqs_category}.")
    
    # exclude singletons if
    if exclude_singletons and only_singletons:
        raise TypeError("'exclude_singletons' and 'only_singletons' can't be provided at the same time!")
    if exclude_singletons:
        mt = io.recalc_info(mt)
        mt = mt.filter_rows(mt.info.AC>1) 
        n_singletons = n_csqs - mt.count()[0]
        print(f"{n_singletons} singletons were excluded!.")
    if only_singletons:
        mt = io.recalc_info(mt)
        mt = mt.filter_rows(mt.info.AC==1) 
        n_singletons = n_csqs - mt.count()[0]
        print(f"{n_singletons} non-singletons were excluded!.")
    # exprot result
    hl.export_vcf(mt, out_prefix + ".vcf.gz") 
    gene_expr = mt.consequence.vep.worst_csq_by_gene_canonical.gene_id
    
    if aggr_method in "fast":
        genes = ko.aggr_phase_count_by_expr(mt, gene_expr)
        genes.checkpoint(out_prefix+"_checkpoint1.mt")
        expr_pko = ko.calc_prob_ko(genes.hom_alt_n, genes.phased, genes.unphased, only_homs=False)
        expr_ko = ko.annotate_knockout(genes.hom_alt_n, expr_pko)
    elif aggr_method in "fast_012":
         genes = ko.aggr_phase_count_by_expr(mt, gene_expr)
         expr_pko = ko.calc_frac_haplotypes(genes.hom_alt_n, genes.phased, genes.unphased)
         expr_ko = ko.annotate_knockout(genes.hom_alt_n, expr_pko)
    elif aggr_method in "only_homs":
         genes = ko.aggr_phase_count_by_expr(mt, gene_expr)
         expr_pko = ko.calc_prob_ko(genes.hom_alt_n, genes.phased, genes.unphased, only_homs=True)
         expr_ko = ko.annotate_knockout(genes.hom_alt_n, expr_pko)
    elif aggr_method in "only_chets":
         genes = ko.aggr_phase_count_by_expr(mt, gene_expr)
         expr_pko = ko.calc_prob_ko(genes.hom_alt_n, genes.phased, genes.unphased, only_chets=True)
         expr_ko = ko.annotate_knockout(genes.hom_alt_n, expr_pko)
    elif aggr_method in "collect":
         genes = ko.collect_phase_count_by_expr(mt, gene_expr)
         genes = ko.sum_gts_entries(genes)
         expr_pko = ko.calc_prob_ko(genes.hom_alt_n, genes.phased, genes.unphased)
         expr_ko = ko.annotate_knockout(genes.hom_alt_n, expr_pko, genes.phased)
    else:
         raise TypeError(str(aggr_method) + " is not allowed!")
    
    # get genen matrix
    genes = genes.annotate_entries(pKO=expr_pko, knockout=expr_ko)

    # setup sites and alleles
    rows = genes.count()[0]
    locus = [ "chr%s:%s" % (chrom, str(i+1)) for i in range(rows)]
    ref = [ get_tid(4) for i in range(rows)]
    alt = [ get_tid(4) for i in range(rows)]

    # convert to HailTable
    df = pd.DataFrame({'locus':locus, 'ref':ref, 'alt':alt})
    ht = hl.Table.from_pandas(df)
    ht = ht.add_index()
    ht = ht.key_by('idx')

    # annotate knockout matrix
    prob = genes.annotate_entries(DS=genes.pKO * 2)
    prob = prob.select_entries(prob.DS)
    prob = prob.add_row_index()
    prob = prob.annotate_rows(
        locus = ht[prob.row_idx].locus,
        tmp_ref = ht[prob.row_idx].ref, 
        tmp_alt = ht[prob.row_idx].alt,
        rsid=prob.gene_id
    )

    # conver to hail locus
    prob = prob.annotate_rows(
            locus=hl.parse_locus(prob.locus),
            alleles=[prob.tmp_ref, prob.tmp_alt]
    )

    # clean up and key appropiately
    prob = prob.key_rows_by(prob.locus, prob.alleles)
    prob = prob.drop(*['gene_id','tmp_ref','tmp_alt','row_idx'])

    # remove invariant sites
    prob = prob.annotate_rows(stdev = hl.agg.stats(prob.DS).stdev)
    prob = prob.filter_rows(prob.stdev > 0)

    # write matrix-table which contains dosages. VCF with only
    # dosages can't be re-read in HAIL, so we write a MatrixTable
    if out_type not in "mt":
        prob = prob.checkpoint(out_prefix + ".mt", overwrite=True)

    # write out variants involved and vcf
    io.export_table(prob, out_prefix, out_type)
    if not only_vcf:
        genes = genes.filter_entries(hl.is_defined(genes.knockout)).entries()
        if aggr_method == "collect":
            genes = genes.transmute(
                        gts=hl.delimit(genes.gts, ";"),
                        varid=hl.delimit(genes.varid, ";")
                        )
        if export_all_gts:
            genes.flatten().export(out_prefix + "_all.tsv.gz")
        else:
            genes = genes.filter(genes.pKO > 0)
            genes.flatten().export(out_prefix + ".tsv.gz")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--chrom', default=None, help='Chromosome to be used') 
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--repartition', default=None, help='Hail repartition files')
    parser.add_argument('--only_vcf', default=False, action='store_true', help='Only return VCF (less memory required when running)')
    parser.add_argument('--checkpoint', default=False, action='store_true', help='Checkpoint gene-aggregation matrix to avoid Spark Memory overflow errors') 
    parser.add_argument('--aggr_method', default=None, required=True, help='How should the CH matrix be generated?')
    # filtering options
    parser.add_argument('--exclude_singletons', default=False, action='store_true', help='Excludes all MAC=1 variants')
    parser.add_argument('--only_singletons', default=False, action='store_true', help='Excludes MAC!=1 variants')
    parser.add_argument('--export_all_gts', default=False, action='store_true', help='Exports a table of all csqs')
    parser.add_argument('--discard_prob_dosages', default=False, action='store_true', help='Discard any dosages < 2 by setting them to zero.')
    parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='What categories should be subsetted to?')

    args = parser.parse_args()

    main(args)



