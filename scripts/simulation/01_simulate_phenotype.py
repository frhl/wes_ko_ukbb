#!/usr/bin/env python3

import hail as hl
import numpy as np
import argparse

from ko_utils import io
from ko_utils import ko
from ukb_utils import hail_init
from ukb_utils import variants
import scipy.stats as stats



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
    prune_hom_alt = args.prune_hom_alt
    seed = args.seed
    h2_nc = try_param_h2(args.h2_nc)
    h2_co = try_param_h2(args.h2_co)
    h2_ko = try_param_h2(args.h2_ko)
    pi_nc = try_param_pi(args.pi_nc)
    pi_co = try_param_pi(args.pi_co)
    pi_ko = try_param_pi(args.pi_ko)
    alpha = try_param_h2(args.alpha)
    beta = try_param_h2(args.beta)
    theta = try_param_h2(args.theta)
    K = float(args.K)
    max_maf = args.max_maf

    # import table
    hail_init.hail_bmrc_init('logs/hail/absence_of_effect.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
    hl.set_global_seed(int(seed))
    mt = io.import_table(in_prefix, in_type)
    
    # filter on MAF 
    if max_maf:
        mt = mt.annotate_rows(MAF=variants.get_maf_expr(mt))
        mt = mt.filter_rows(mt.MAF <= float(max_maf))

    # annotate with variant consequence
    mt = mt.explode_rows(mt.consequence.vep.worst_csq_by_gene_canonical)
    mt = mt.annotate_rows(
        consequence_category=ko.csqs_case_builder(
                worst_csq_expr=mt.consequence.vep.worst_csq_by_gene_canonical,
                use_loftee=True))

    # filter to variants with SD > 0
    mt = mt.filter_rows(hl.agg.stats(mt.GT.n_alt_alleles()).stdev>0)
    
    # subset to potential 'damaging variants' 
    items_co = ['pLoF','LC','damaging_missense']
    items_nc = ['non_coding'] 
    mt_co = mt.filter_rows(hl.literal(set(items_co)).contains(mt.consequence_category))
    mt_nc = mt.filter_rows(hl.literal(set(items_nc)).contains(mt.consequence_category))

    c1 = mt.count()
    c_nc = mt_nc.count()
    c_co = mt_co.count()
    print("c1 = %s .. nc = %s .. co = %s" % (c1, c_nc, c_co))

    if h2_nc > 0 or alpha:
        # simulate alphas (coding region)
        print("Simulating non-coding variant phenotype with h2=%s and pi=%s" % (h2_nc, pi_nc))
        if alpha:
            alpha = ko.make_effect_size(mt_nc, alpha, pi_nc)
        else:
            alpha = ko.simulate_effect_size(mt_nc, h2_nc, pi_nc)
        mt_nc = mt_nc.annotate_rows(alpha = alpha)
        mt_nc = hl.experimental.ldscsim.normalize_genotypes(mt_nc.GT)
        mt_nc = mt.annotate_cols(y_no_noise_nc=hl.agg.sum(mt_nc.alpha * mt_nc['norm_gt']))
        mt = mt.annotate_cols(y_no_noise_nc = mt_nc.index_cols(mt.col_key).y_no_noise_nc) 

    if h2_co > 0 or beta:
        # simulate betas (coding region)
        print("Simulating coding variant phenotype with h2=%s and pi=%s" % (h2_co, pi_co))
        if beta:
            beta = ko.make_effect_size(mt_co, beta, pi_co)
        else:
            beta = ko.simulate_effect_size(mt_co, h2_co, pi_co)
        mt_co = mt_co.annotate_rows(beta = beta)
        mt_co = hl.experimental.ldscsim.normalize_genotypes(mt_co.GT)
        mt_co = mt_co.annotate_cols(y_no_noise_co=hl.agg.sum(mt_co.beta * mt_co['norm_gt']))
        mt = mt.annotate_cols(y_no_noise_nc = mt_co.index_cols(mt.col_key).y_no_noise_co) 

    # annotate with CH effect
    if h2_ko > 0 or theta:

        # simulate compound het effects
        print("Simulating compound het phenotype with h2=%s and pi=%s" % (h2_ko, pi_ko))

        # prune away knockedout owed to homozygote alternates
        gene_mt = mt_nc
        if prune_hom_alt:
            prune_hom_alt = float(prune_hom_alt)
            gene_mt = gene_mt.transmute_entries(GT = ko.rand_hom_to_het(gene_mt.GT, prune_hom_alt))

        # generate gene x sample matrix
        gene_expr = gene_mt.consequence.vep.worst_csq_by_gene_canonical.gene_id
        genes = ko.aggr_phase_count_by_expr(gene_mt, gene_expr)
        expr_pko = ko.calc_prob_ko(genes.hom_alt_n, genes.phased, genes.unphased)
        expr_ko = ko.annotate_knockout(genes.hom_alt_n, expr_pko)
        genes = genes.annotate_entries(pKO=expr_pko, knockout=expr_ko)

        n = genes.count()[0]
        if pi_ko:
            n_causal_genes = n * pi_ko
            print("Around %0.2f of %0.0f genes will be causal with pi=%f." % (n_causal_genes, n, pi_ko))

        # simulate thetas (CH effects).
        genes = genes.filter_rows(hl.agg.stats(genes.pKO).stdev>0)
        if theta:
            theta = ko.make_effect_size(genes, theta, pi_ko)
        else:
            theta = ko.simulate_effect_size(genes, h2_ko, pi_ko)
        genes = genes.annotate_rows(theta = theta)
        #genes = ko.make_thetas(genes, h2=h2_ko, pi=pi_ko)
        genes = ko.normalize_by_name(genes, "pKO")
        genes = genes.filter_rows(genes.theta > 0)
        genes = genes.annotate_cols(y_no_noise_ko=hl.agg.sum(genes.theta * genes.norm_pKO))
        genes = genes.annotate_cols(y = genes.y_no_noise_ko + hl.rand_norm(0, hl.sqrt(1-h2_ko))) 
    
        if pi_ko:
            found_causal_genes = np.sum(np.array(genes.theta.collect()) != 0)
            print("Observed %0.2f of %0.0f causal genes.." % (found_causal_genes, n))

        # return thetas for genes and samples
        ht = genes.select_rows('theta').select_entries(*['pKO','norm_pKO','knockout'])
        ht.entries().flatten().export(out_prefix + "_genes.tsv.gz")

        # annotate original matrix with thetas from gene x sample matrix
        mt = mt.annotate_cols(y_no_noise_ko = genes.index_cols(mt.col_key).y_no_noise_ko) 


    # generate final phenotype by adding noise component
    if h2_ko == 0 and h2_co == 0 and h2_nc == 0:
           mt = mt.annotate_cols(y_cts=hl.rand_norm(0, hl.sqrt(1)))
    elif h2_ko > 0 and h2_co > 0 and h2_nc == 0:
        mt = mt.annotate_cols(y_cts=mt.y_no_noise_co +
                    mt.y_no_noise_ko + hl.rand_norm(0, hl.sqrt(1-h2_ko-h2_co)))
    elif h2_ko > 0 and h2_co == 0 and h2_nc == 0:
        mt = mt.annotate_cols(y_cts = mt.y_no_noise_ko +
                    hl.rand_norm(0, hl.sqrt(1-h2_ko)))
    elif h2_ko == 0 and h2_co > 0 and h2_nc == 0:
        mt = mt.annotate_cols(y_cts = mt.y_no_noise_co +
                    hl.rand_norm(0, hl.sqrt(1-h2_co)))
    elif h2_ko > 0 and h2_co > 0 and h2_nc > 0:
        mt = mt.annotate_cols(y_cts=mt.y_no_noise_co +
                    mt.y_no_noise_ko + hl.rand_norm(0, hl.sqrt(1-h2_ko-h2_co-h2_nc)))
    elif h2_ko > 0 and h2_co == 0 and h2_nc > 0:
        mt = mt.annotate_cols(y_cts = mt.y_no_noise_ko +
                    hl.rand_norm(0, hl.sqrt(1-h2_ko-h2_nc)))
    elif h2_ko == 0 and h2_co > 0 and h2_nc > 0:
        mt = mt.annotate_cols(y_cts = mt.y_no_noise_co +
                    hl.rand_norm(0, hl.sqrt(1-h2_co-h2_nc)))
    else:
        raise TypeError("Invalid use of h2_ko, h2_co and h2_nc! Are some of them None?")

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
    parser.add_argument('--causal_genes_ko', default=None, help='Desired number of causal genes')
    parser.add_argument('--K', default=0, help='Prevalence of phenotype: cases / (cases + controls)')
    parser.add_argument('--max_maf', default=None, help='Maximum minor allele frequency')
    parser.add_argument('--seed', default=None, help='seed for random simulations')
    parser.add_argument('--in_prefix', default=None, help='Path prefix for input dataset')
    parser.add_argument('--in_type', default=None, help='Either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--prune_hom_alt', default=None, help='Path prefix for output dataset')
    parser.add_argument('--alpha', default=None, help='Pre-defined effect size for causal non-coding variants')
    parser.add_argument('--beta', default=None, help='Pre-defined effect size for causal coding variants')
    parser.add_argument('--theta', default=None, help='Pre-defined effect size for causal compound het pseudo variants')
    #parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='comma sepearted field')
    #parser.add_argument('--out_type', default=None, help='Either "mt", "vcf" or "plink"')
    args = parser.parse_args()

    main(args)


