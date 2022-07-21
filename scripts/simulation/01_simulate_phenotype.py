#!/usr/bin/env python3

import hail as hl
import numpy as np
import argparse

from ko_utils import io
from ko_utils import ko
from ukb_utils import hail_init
from ukb_utils import variants
import scipy.stats as stats


def rescale_variances_hail(X, mt, target_variance):
    stats = mt.aggregate_cols(hl.agg.stats(X))
    return(rescale_variances(X, stats.stdev ** 2, target_variance))
    
def rescale_variances(X, X_variance, Y_variance):
    X_sd = X_variance ** (1/2)
    Y_sd = Y_variance ** (1/2)
    K = Y_sd / X_sd
    Y = K*X
    return(Y)

def make_effect(mt, h2, pi = None):
    """ Make effect sizes for either infintesimal or spike and slab model. """
    M = mt.count()[0]
    pi_temp = 1 if pi == None else pi
    return(hl.rand_bool(pi_temp)*hl.rand_norm(0, hl.sqrt(h2/(M*pi_temp))))



def try_param_h2(x):
    return float(x) if x not in [None, "NA","None"] else None

def try_param_pi(x):
    return float(x) if x not in [0, None, "NA","None", "0.0", "0.00", "0.000", "0"] else None

class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))

def main(args):

    in_prefix = args.in_prefix
    in_type = args.in_type
    out_prefix = args.out_prefix
    seed = args.seed
    max_maf = args.max_maf
    h2 = try_param_h2(args.h2)
    var_beta = try_param_h2(args.var_beta)
    var_theta = try_param_h2(args.var_theta)
    pi_beta = try_param_pi(args.pi_beta)
    pi_theta = try_param_pi(args.pi_theta)
    K = float(args.K)

    # import table
    hail_init.hail_bmrc_init('logs/hail/simulate_phenotype.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
    hl.set_global_seed(int(seed))
    mt = io.import_table(in_prefix, in_type)
    
    # filter on MAF 
    if max_maf:
        mt = mt.annotate_rows(MAF=variants.get_maf_expr(mt))
        mt = mt.filter_rows(mt.MAF <= float(max_maf))

    # annotate with variant consequence and collapse by gene
    mt = mt.explode_rows(mt.consequence.vep.worst_csq_by_gene_canonical)
    mt = mt.annotate_rows(
        consequence_category=ko.csqs_case_builder(
                worst_csq_expr=mt.consequence.vep.worst_csq_by_gene_canonical,
                use_loftee=True))

    # filter to variants with SD > 0 and codign variants
    categories = ['pLoF','LC','damaging_missense']
    mt = mt.filter_rows(hl.agg.stats(mt.GT.n_alt_alleles()).stdev>0)
    mt = mt.filter_rows(hl.literal(set(categories)).contains(mt.consequence_category))

    # collapse to gene level
    gene_expr = mt.consequence.vep.worst_csq_by_gene_canonical.gene_id
    mt = ko.aggr_phase_count_by_expr(mt, gene_expr)

    # annotate knockout type
    expr_pko = ko.calc_prob_ko(mt.hom_alt_n, mt.phased, mt.unphased)
    expr_ko = ko.annotate_knockout(mt.hom_alt_n, expr_pko)
    mt = mt.annotate_entries(
        pKO=expr_pko, 
        knockout=expr_ko
    )   

    # What haplotypes are affected?
    mt = mt.annotate_entries(
        H1 = 1*((mt.phased.a1 > 0) | (mt.hom_alt_n > 0)),
        H2 = 1*((mt.phased.a2 > 0) | (mt.hom_alt_n > 0)),
    )

    # combine into single dosage (G) matrix 
    mt = mt.annotate_entries(G=hl.int32(mt.H1+mt.H2))
    mt = mt.annotate_rows(**{'stats': hl.agg.stats(mt.G)})
    mt = mt.filter_rows(mt.stats.stdev > 0)
    mt = mt.annotate_entries(
        G_norm=(mt.G-mt.stats.mean)/mt.stats.stdev
    )
    
    print("h2:", str(h2))
    print("var_beta:", str(var_beta))
    print("var_theta:", str(var_theta))
    # setup effects
    mt = mt.annotate_rows(beta = make_effect(mt, h2=var_beta, pi=pi_beta))
    mt = mt.annotate_rows(theta_nosign = make_effect(mt, h2=var_theta, pi=pi_theta))

    # keep the sign to ensure that effects are always in the same direction
    mt = mt.annotate_rows(
        beta_sign = hl.case().when(hl.sign(mt.beta) == 0, 1).default(hl.sign(mt.beta))
    )

    # we would like to operate on the same genotype scale as the additive effects,
    # so use the original scaled matrix, but only keep hom/ch alt alleles sites
    mt = mt.annotate_entries(
        G_norm_alt = (hl.case().when(mt.G == 2, mt.G_norm).default(0))
    )

    # set sign
    mt = mt.annotate_rows(
        theta = mt.theta_nosign * mt.beta_sign
    )

    # sum up effect from domincance
    mt = mt.annotate_cols(y_no_noise_add=hl.agg.sum(mt.beta * mt.G_norm))
    mt = mt.annotate_cols(y_no_noise_dom=hl.agg.sum(mt.theta * mt.G_norm_alt))
    mt = mt.annotate_cols(y_no_noise=mt.y_no_noise_add+mt.y_no_noise_dom)

    # re-scale effects
    if h2 > 0 and (var_beta + var_theta) > 0:
        mt = mt.annotate_cols(y_no_noise_rescaled = rescale_variances_hail(mt.y_no_noise, mt, h2))
    else:
        mt = mt.annotate_cols(y_no_noise_rescaled = 0)    

    # add up all effects
    mt = mt.annotate_cols(y_noise = hl.rand_norm(0, hl.sqrt(1-h2)))
    mt = mt.annotate_cols(y = mt.y_noise + mt.y_no_noise_rescaled)
    mt = mt.annotate_cols(dom_add_frac = mt.y_no_noise_dom / (mt.y_no_noise_dom + mt.y_no_noise_add))

    # binarize phenotype
    if K is not None:
        y_stats = mt.aggregate_cols(hl.agg.stats(mt.y))
        threshold = stats.norm.ppf(1-K, loc=y_stats.mean, scale=y_stats.stdev)
        mt = mt.annotate_cols(case=mt.y > threshold)

    # export effect sizes
    ht = mt.select_rows(*['beta','theta']).select_entries(*['pKO','knockout'])
    ht.entries().flatten().export(out_prefix + "_entries.tsv.gz")

    # export simulated phenotypes
    ht = mt.cols()
    ht.flatten().export(out_prefix + "_cols.tsv.gz")

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('--chrom', default=None, help='chromosome')
    parser.add_argument('--h2', default=0, help='Heritability for coding region')
    parser.add_argument('--var_beta', default=0, help='Heritability for coding region')
    parser.add_argument('--var_theta', default=0, help='Heritability for knockouts')
    parser.add_argument('--pi_beta', default=0, help='Probability of variant being causal in coding region')
    parser.add_argument('--pi_theta', default=0, help='Probability of variant being causal in non-coding region')
    parser.add_argument('--K', default=0, help='Prevalence of phenotype: cases / (cases + controls)')
    parser.add_argument('--max_maf', default=None, help='Maximum minor allele frequency')
    parser.add_argument('--seed', default=None, help='seed for random simulations')
    parser.add_argument('--in_prefix', default=None, help='Path prefix for input dataset')
    parser.add_argument('--in_type', default=None, help='Either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


