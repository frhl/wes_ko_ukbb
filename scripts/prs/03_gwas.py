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
   
    if phenotypes:
        ht = hl.import_table(phenotypes,
                     types={'eid': hl.tstr},
                     missing=["",'""',"NA"],
                     impute=True,
                     force=True,
                     key='eid'
                     ) 
        mt = mt.annotate_cols(pheno=ht[mt.s])
    else:
        raise ValueError("param 'phenotypes' is not set! ")

    split_covariates = covariates.split(',')
    if set(list(mt.pheno)) >= set(split_covariates):
        covariates = [mt.pheno[x] for x in split_covariates]
        covariates.insert(0, 1) 
        if response in list(mt.pheno):
            if min_maf_cutoff:
               mt = mt.filter_rows(variants.get_maf_expr(mt) > float(min_maf_cutoff))
            if mt.pheno[response].dtype == hl.dtype('float64'):
               reg = hl.linear_regression_rows(
                        y=mt.pheno[response],
                        x=mt.GT.n_alt_alleles(),
                        covariates=covariates
                        )
            elif mt.pheno[response].dtype == hl.dtype('bool'):
                if response not in list(mt.pheno):
                    raise ValueError("Response: " + str(response) + " is not in pheno file!")
                n_total = mt.aggregate_cols(hl.agg.sum(hl.is_defined(mt.pheno[response])))
                n_cases = mt.aggregate_cols(hl.agg.sum(mt.pheno[response] == True))
                n_controls = n_total - n_cases
                #controls = mt.aggregate_cols(hl.agg.sum(mt.pheno[response] == False))
                if n_total < 1:
                    raise ValueError(str(defined) + " cases/controls defined for phenotype " + str(response)+". Exiting..")
                if n_cases < int(min_cases):
                     raise ValueError(str(n_cases) + " cases found! Expected +" + str(min_cases))
                if adjust_maf_by_case_control:
                     min_maf = hl.max(0.01, 25/(2 * hl.min([n_cases, n_controls]))).collect()[0]
                     mt = mt.filter_rows(variants.get_maf_expr(mt) > min_maf)
                     # need to re-initalise covariates, otherwise
                     # hail will throw an error about using multiple mts
                     covariates = [mt.pheno[x] for x in split_covariates]
                     covariates.insert(0, 1) 
                     with open("data/prs/sumstat/binary/maf_thresholds.txt","a") as outfile:
                         outfile.write(f"{input_path}\t{response}\t{n_cases}\t{n_controls}\t{min_maf}\n")
                
                reg = hl.logistic_regression_rows(
                        test='wald',
                        y=mt.pheno[response],
                        x=mt.GT.n_alt_alleles(),
                        covariates=covariates
                        )
                # get effective sample size
                #n_total = mt.aggregate_cols(hl.agg.sum(hl.is_defined(mt.pheno[response])))
                #n_cases = mt.aggregate_cols(hl.agg.sum(mt.pheno[response] == True))
                reg = reg.annotate(n=n_total, n_cases=n_cases, maf_cutoff = min_maf)
            else:
                raise TypeError("Response variable is not a float64 or boolean!")
        else:
            raise TypeError("Response variable is not in phenotype file")
    else:    
        raise TypeError("Some or all covariates are not in phenotype file!")

    # get AFs
    reg = reg.annotate(
        AN=mt.rows()[reg.key].info.AN,
        AC=mt.rows()[reg.key].info.AC,
        AF=mt.rows()[reg.key].info.AF
    )
    
    # Get individual alleles
    reg = reg.annotate(
        a0=reg.alleles[0],
        a1=reg.alleles[1]
    )

    reg.flatten().export(out_prefix + ".txt.gz")



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


