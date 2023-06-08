

#!/usr/bin/env python3

import hail as hl
import numpy as np
import argparse

from ko_utils import io
from ko_utils import ko
from ukb_utils import hail_init
from ukb_utils import variants
import scipy.stats as stats


def make_effect(mt, h2, b, pi = None):
    """ Make effect sizes for either infintesimal or spike and slab model. """
    M = mt.count()[0]
    pi_temp = 1 if pi == None else pi
    mean = b/(M*pi) 
    variance = h2/(M*pi)
    print("h2=" + str(h2) + " b=" + str(b) + " mean=" + str(mean) + "variance=" + str(variance))
    return(hl.rand_bool(pi_temp)*hl.rand_norm(mean, hl.sqrt(variance)))

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
    h2 = float(args.h2) #try_param_h2(args.h2)
    b = float(args.b) #try_param_h2(args.b)
    pi = float(args.pi)
    K = float(args.K)

    # import table
    hail_init.hail_bmrc_init('logs/hail/simulate_phenotype.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
    hl.set_global_seed(int(seed))
    mt = io.import_table(in_prefix, in_type)
    
    # setup effects
    if h2 > 0:
        # both G_add_norm and G_rec_norm are on the same scale, however G_rec_norm
        # has been normalized AFTER excluding DS < 2, i.e. only recessive effects will be allowed.
        mt = mt.annotate_rows(theta=make_effect(mt, h2=h2, b=b, pi=pi))
        mt = mt.annotate_cols(y_no_noise=hl.agg.sum(mt.theta * mt.G_rec_norm_by_add))
    else:
        mt = mt.annotate_cols(y_no_noise = 0)
        mt = mt.annotate_rows(
                theta = 0,
       )

    # add up all effects
    mt = mt.annotate_cols(y_noise = hl.rand_norm(0, hl.sqrt(1-h2)))
    mt = mt.annotate_cols(y = mt.y_noise + mt.y_no_noise)

    # binarize phenotype
    if K is not None:
        y_stats = mt.aggregate_cols(hl.agg.stats(mt.y))
        threshold = stats.norm.ppf(1-K, loc=y_stats.mean, scale=y_stats.stdev)
        mt = mt.annotate_cols(case=mt.y > threshold)

    # export effect sizes
    ht = mt.select_rows(*['rsid','theta']).select_entries(*['pKO','knockout'])
    ht.entries().flatten().export(out_prefix + "_entries.tsv.gz")

    # export simulated phenotypes
    ht = mt.cols()
    ht.flatten().export(out_prefix + "_cols.tsv.gz")

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('--chrom', default=None, help='chromosome')
    parser.add_argument('--h2', default=0, help='Heritability for trait')
    parser.add_argument('--b', default=0, help='Average effect for a variant')
    parser.add_argument('--pi', default=0, help='Probability of variant being causal')
    parser.add_argument('--K', default=0, help='Prevalence of phenotype: cases / (cases + controls)')
    parser.add_argument('--seed', default=None, help='seed for random simulations')
    parser.add_argument('--in_prefix', default=None, help='Path prefix for input dataset')
    parser.add_argument('--in_type', default=None, help='Either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


