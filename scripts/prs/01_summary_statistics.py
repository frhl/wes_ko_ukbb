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
    phenotypes = args.phenotypes
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    
    if dataset in "imp":
        mt = genotypes.get_ukb_imputed_v3_bgen(chroms=[chrom])
        if min_info:
            ht = genotypes.get_ukb_parsed_imputed_v3_mfi(chroms=[chrom])
            mt = mt.annotate_rows(info_score = ht[(mt.locus, mt.alleles)].info)
            mt = mt.filter_rows(mt.info_score >= float(min_info))
    elif dataset in "calls":
        mt = genotypes.get_ukb_genotypes_bed(chroms=[chrom])
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
                     impute = True, 
                     key = 'eid') 
        mt = mt.annotate_cols(pheno = ht[mt.s])
    else:
        raise ValueError("param 'phenotypes' is not set! ")
  

    # Read file of covariates 

    # check that covariates are actually in list(mt.pheno)

    # response variable

    # create list of explanatory variables
    explanatory = [1, mt.GT.n_alt_alleles(), ... ]

    current_pheno = list(mt.pheno)[index]
    
    mt = mt.annotate_rows(linreg=hl.agg.linreg(y=mt.pheno[current_pheno],x=[1,mt.GT.n_alt_alleles()]))
    
    #hl.linear_regression_rows(y=mt.pheno[''],
    #    x=mt.GT.n_alt_alleles(),
    #    covariates=[1])

    #qc.export_table(mt, out_prefix, out_type)


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='chromosome')
    parser.add_argument('--min_info', default=None, help='minimum info score (only imputed data)')
    parser.add_argument('--liftover', default=None, action='store_true', help='perform liftover')
    parser.add_argument('--extract_samples', default=None, help='Subset to sample IDs in file')
    parser.add_argument('--dataset', default=None, help='Either "imp" or "calls".')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Either "mt", "vcf" or "plink"')
    args = parser.parse_args()

    main(args)


