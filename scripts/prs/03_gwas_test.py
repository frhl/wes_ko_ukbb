#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import variants

def main(args):

    # parser
    input_path = args.input_path
    input_type = args.input_type
    adjust_maf_by_case_control = args.adjust_maf_by_case_control
    response = args.response
    min_cases = args.min_cases
    min_maf_cutoff = args.min_maf_cutoff
    covariates = args.covariates
    phenotypes = args.phenotypes
    out_prefix = args.out_prefix

    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    mt = io.import_table(input_path, input_type)

    hl.init(log='hail.log')
    mt = hl.balding_nichols_model(1,100,2)
    mt = mt.annotate_cols(y = hl.rand_bool(0.5))
    result_ht = hl.logistic_regression_rows(
        test='wald',
        y=mt.y, 
        x=mt.GT.n_alt_alleles(), 
        covariates=[1]
    )   

    result_ht.show()

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='chromosomes to be parsed.')
    parser.add_argument('--input_path', default=None, help='Path to input.')
    parser.add_argument('--input_type', default=None, help='Input type (vcf/mt/plink).')
    parser.add_argument('--min_cases', default=1, help='Minimum number of cases allowed for binary traits.')
    parser.add_argument('--min_maf_cutoff', default=None, help='filter by minimum minor allele frequency')
    parser.add_argument('--adjust_maf_by_case_control', default=None, action='store_true', help='Perform further subsets on MAF based on case/controls')
    #parser.add_argument('--liftover', default=None, action='store_true', help='perform liftover')
    #parser.add_argument('--extract_samples', default=None, help='Subset to sample IDs in file')
    #parser.add_argument('--hapmap', default=None, help='Path to HapMap SNPs')
    #parser.add_argument('--dataset', default=None, help='Either "imp" or "calls".')
    parser.add_argument('--response', default=None, help='Response variable')
    parser.add_argument('--covariates', default=None, help='list of covariates')
    parser.add_argument('--phenotypes', default=None, help='File path to phenotypes.')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


