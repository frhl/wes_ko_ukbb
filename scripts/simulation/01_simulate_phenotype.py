#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ko_utils import ko
from ukb_utils import hail_init
import scipy.stats as stats


def try_set_float(x):
        return float(x) if x not in [None, "NA","None"] else None

def make_y_region(mt: hl.MatrixTable, h2, pi, region: list()):
    mt = mt.filter_rows(hl.literal(set(region)).contains(mt.consequence_category))
    mt = mt.filter_rows(hl.agg.stats(mt.GT.n_alt_alleles()).stdev>0)
    mt = hl.experimental.ldscsim.make_betas(mt, h2=h2, pi=pi)[0]
    mt = hl.experimental.ldscsim.normalize_genotypes(mt.GT) 
    mt = mt.annotate_cols(y_no_noise=hl.agg.sum(mt.beta * mt['norm_gt']))
    return(mt)

class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))

def main(args):

    in_prefix = args.in_prefix
    in_type = args.in_type
    out_prefix = args.out_prefix
    prune_hom_alt = args.prune_hom_alt
    seed = args.seed
    h2_nc = try_set_float(args.h2_nc)
    h2_co = try_set_float(args.h2_co)
    h2_ko = try_set_float(args.h2_ko)
    pi_nc = try_set_float(args.pi_nc)
    pi_co = try_set_float(args.pi_co)
    pi_ko = try_set_float(args.pi_ko)
    K = try_set_float(args.K)

    # import table
    hail_init.hail_bmrc_init('logs/hail/absence_of_effect.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
    hl.set_global_seed(int(seed))
    mt = io.import_table(in_prefix, in_type)

    # annotate with variant consequence
    mt = mt.explode_rows(mt.consequence.vep.worst_csq_by_gene_canonical)
    mt = mt.annotate_rows(
        consequence_category=ko.csqs_case_builder(
                worst_csq_expr=mt.consequence.vep.worst_csq_by_gene_canonical,
                use_loftee=True))

    # subset to potential 'damaging variants' 
    items = ['pLoF','LC','damaging_missense']
    mt = mt.filter_rows(hl.literal(set(items)).contains(mt.consequence_category))
 
   

    # simulate betas (coding region)
    mt = mt.filter_rows(hl.agg.stats(mt.GT.n_alt_alleles()).stdev>0)
    mt = hl.experimental.ldscsim.make_betas(mt, h2=h2_snp, pi=pi_snp)[0]
    mt = hl.experimental.ldscsim.normalize_genotypes(mt.GT) 
    mt = mt.annotate_cols(y_no_noise_snp=hl.agg.sum(mt.beta * mt['norm_gt']))

    # annotate with CH effect
    if h2_ko is not None:

        gene_mt = mt

        # prune away knockedout owed to homozygote alternates
        if prune_hom_alt:
            prune_hom_alt = float(prune_hom_alt)
            gene_mt = gene_mt.transmute_entries(GT = ko.rand_hom_to_het(gene_mt.GT, prune_hom_alt))

        # generate gene x sample matrix
        gene_expr = gene_mt.consequence.vep.worst_csq_by_gene_canonical.gene_id
        genes = ko.aggr_phase_count_by_expr(gene_mt, gene_expr)
        expr_pko = ko.calc_prob_ko(genes.hom_alt_n, genes.phased, genes.unphased)
        expr_ko = ko.annotate_knockout(genes.hom_alt_n, expr_pko)
        genes = genes.annotate_entries(pKO=expr_pko, knockout=expr_ko)

        # simulate thetas (CH effects).
        genes = genes.filter_rows(hl.agg.stats(genes.pKO).stdev>0)
        genes = ko.normalize_by_name(genes, "pKO") # <---- set library
        genes = ko.make_thetas(genes, h2=h2_ko, pi=pi_ko)
        genes = genes.annotate_cols(y_no_noise_ko=hl.agg.sum(genes.theta * genes.norm_pKO))
        
        # return thetas for genes and samples
        ht = genes.select_rows('theta').select_entries(*['pKO','knockout'])
        ht.entries().flatten().export(out_prefix + "_genes.tsv.gz")

        # annotate original matrix with thetas from gene x sample matrix
        mt = mt.annotate_cols(y_no_noise_ko = genes.index_cols(mt.col_key).y_no_noise_ko)
        mt = mt.annotate_cols(y_cts=mt.y_no_noise_snp +
                mt.y_no_noise_ko + hl.rand_norm(0, hl.sqrt(1-h2_ko-h2_snp)))
    else:
        mt = mt.annotate_cols(y_cts = mt.y_no_noise_snp +
                hl.rand_norm(0, hl.sqrt(1-h2_snp)))

    # binarize phenotype
    if K is not None:
        y_stats = mt.aggregate_cols(hl.agg.stats(mt.y_cts))
        threshold = stats.norm.ppf(1-K, loc=y_stats.mean, scale=y_stats.stdev)
        mt = mt.annotate_cols(y_bin=mt.y_cts > threshold) 

    # export simulated phenotypes
    ht = mt.cols()
    ht.flatten().export(out_prefix + ".tsv.gz")

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('--chrom', default=None, help='chromosome')
    parser.add_argument('--h2_co', default=0, help='Heritability for coding region')
    parser.add_argument('--h2_nc', default=0, help='Heritability for non-coding region')
    parser.add_argument('--h2_ko', default=0, help='Heritability for knockouts')
    parser.add_argument('--pi_co', default=0, help='Probability of variant being causal in coding region')
    parser.add_argument('--pi_nc', default=0, help='Probability of variant being causal in non-coding region')
    parser.add_argument('--pi_ko', default=0, help='Probability of variant being causal in knockouts')
    parser.add_argument('--K', default=0, help='Prevalence of phenotype: cases / (cases + controls)')
    parser.add_argument('--seed', default=None, help='seed for random simulations')
    parser.add_argument('--in_prefix', default=None, help='Path prefix for input dataset')
    parser.add_argument('--in_type', default=None, help='Either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--prune_hom_alt', default=None, help='Path prefix for output dataset')
    #parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='comma sepearted field')
    #parser.add_argument('--out_type', default=None, help='Either "mt", "vcf" or "plink"')
    args = parser.parse_args()

    main(args)


