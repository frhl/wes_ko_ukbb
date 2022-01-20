#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import qc
from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import variants
from ukb_utils import samples


def main(args):

    # parser
    chrom = args.chrom
    dataset = args.dataset
    extract_samples = args.extract_samples
    min_info = args.min_info
    liftover = args.liftover
    response = args.response
    covariates = args.covariates
    phenotypes = args.phenotypes
    out_prefix = args.out_prefix

    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    
    if dataset in "imp":
        mt = genotypes.get_ukb_imputed_v3_bgen(chroms=chrom)
        if min_info:
            ht = genotypes.get_ukb_parsed_imputed_v3_mfi(chroms=chrom)
            mt = mt.annotate_rows(info_score = ht[(mt.locus, mt.alleles)].info)
            mt = mt.filter_rows(mt.info_score >= float(min_info))
    elif dataset in "calls":
        mt = genotypes.get_ukb_genotypes_bed(chroms=chrom)
    else:
        raise TypeError(f"{dataset} is not 'imp' or 'calls'")
    
    if extract_samples:
        ht_samples = hl.import_table(extract_samples, no_header=True, key='f0', delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_samples[mt.col_key]))
    
    if liftover:
        mt = variants.liftover(mt, from_build='GRCh37', to_build='GRCh38', drop_annotations=True)
   
    if phenotypes:
        ht = hl.import_table(phenotypes,
                     types={'eid': hl.tstr},
                     missing=["",'""'],
                     impute=True,
                     key='eid') 
        mt = mt.annotate_cols(pheno=ht[mt.s])
    else:
        raise ValueError("param 'phenotypes' is not set! ")
    
    covariates = covariates.split(',')
    if set(list(mt.pheno)) >= set(covariates):
        covariates.insert(0, 1)
        if response in list(mt.pheno):
            if mt.pheno[response].dtype == hl.dtype('float64'):
                reg = hl.linear_regression_rows(
                        x=mt.pheno[response],
                        y=mt.GT.n_alt_alleles(),
                        covariates=covariates
                        )
            elif mt.pheno[response].dtype == hl.dtype('bool'):
                reg = hl.logistic_regression_rows(
                        test='wald',
                        y=mt.pheno[response],
                        x=mt.GT.n_alt_alleles(),
                        covariates=covariates
                        )
            else:
                raise TypeError("Response variable is not a float64 or boolean!")
        else:
            raise TypeError("Response variable is not in phenotype file")
    else:    
        raise TypeError("Some or all covariates are not in phenotype file!")

    reg.flatten.export(out_prefix + ".txt.gz")



if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, nargs='+', help='chromosomes to be parsed.')
    parser.add_argument('--min_info', default=None, help='minimum info score (only imputed data)')
    parser.add_argument('--liftover', default=None, action='store_true', help='perform liftover')
    parser.add_argument('--extract_samples', default=None, help='Subset to sample IDs in file')
    parser.add_argument('--dataset', default=None, help='Either "imp" or "calls".')
    parser.add_argument('--response', default=None, help='Response variable')
    parser.add_argument('--covariates', default=None, help='list of covariates')
    parser.add_argument('--phenotypes', default=None, help='File path to phenotypes.')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


