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
    """ Rescale variance X in of field in mt by some target variance
    """
    stats = mt.aggregate_cols(hl.agg.stats(X))
    return(rescale_variances(X, stats.stdev ** 2, target_variance))
    
def rescale_variances(X, X_variance, Y_variance):
    """ Variance rescaling equation
    """
    X_sd = X_variance ** (1/2)
    Y_sd = Y_variance ** (1/2)
    K = Y_sd / X_sd
    Y = K*X
    return(Y)

def make_effect(mt, h2, pi = None, causal_bool = None):
    """ Make effect sizes for either infintesimal or spike and slab model. """
    M = mt.count()[0]
    pi_temp = 1 if pi == None else pi
    if causal_bool is not None:
        return(causal_bool*hl.rand_norm(0, hl.sqrt(h2/(M*pi_temp))))
    else:
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
    
    print("h2:", str(h2))
    print("var_beta:", str(var_beta))
    print("var_theta:", str(var_theta))
    
    # setup effects
    if h2 > 0 and (var_beta + var_theta) > 0:

        causal_pi = 0.2
        causal_bool = hl.rand_bool(causal_pi)

        # setup standard additive effects
        mt = mt.annotate_rows(beta = make_effect(mt, h2=var_beta, pi=causal_pi, causal_bool=causal_bool))

        # keep track of sign for additive effects
        mt = mt.annotate_rows(
            beta_sign = hl.case().when(hl.sign(mt.beta) == 0, 1).default(hl.sign(mt.beta))
        )

        # setup recessive effects
        mt = mt.annotate_rows(
                theta_nosign = hl.abs(make_effect(mt, h2=var_theta, pi=causal_pi, causal_bool=causal_bool))
        )
        
        # keep theta sign consistent with beta sign
        mt = mt.annotate_rows(
            theta = mt.theta_nosign * mt.beta_sign
        )

        # both G_add_norm and G_rec_norm are on the same scale, however G_rec_norm
        # has been normalized AFTER excluding DS < 2, i.e. only recessive effects will be allowed.
        mt = mt.annotate_cols(y_no_noise_add=hl.agg.sum(mt.beta * mt.G_add_norm)) # beta
        mt = mt.annotate_cols(y_no_noise_rec=hl.agg.sum(mt.theta * mt.G_rec_norm)) # theta
        mt = mt.annotate_cols(y_no_noise=mt.y_no_noise_add+mt.y_no_noise_rec)

        # re-scale effects genetic effects accordingly 
        mt = mt.annotate_cols(y_no_noise_rescaled = rescale_variances_hail(mt.y_no_noise, mt, h2))
    else:
        mt = mt.annotate_cols(y_no_noise_rescaled = 0)
        mt = mt.annotate_rows(
                beta = 0,
                theta = 0
       )

    # add up all effects
    mt = mt.annotate_cols(y_noise = hl.rand_norm(0, hl.sqrt(1-h2)))
    mt = mt.annotate_cols(y = mt.y_noise + mt.y_no_noise_rescaled)

    # binarize phenotype
    if K is not None:
        y_stats = mt.aggregate_cols(hl.agg.stats(mt.y))
        threshold = stats.norm.ppf(1-K, loc=y_stats.mean, scale=y_stats.stdev)
        mt = mt.annotate_cols(case=mt.y > threshold)

    # export effect sizes
    ht = mt.select_rows(*['rsid','beta','theta']).select_entries(*['pKO','knockout'])
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
    parser.add_argument('--seed', default=None, help='seed for random simulations')
    parser.add_argument('--in_prefix', default=None, help='Path prefix for input dataset')
    parser.add_argument('--in_type', default=None, help='Either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


