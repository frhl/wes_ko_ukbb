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
    min_maf_cutoff = args.min_maf_cutoff
    covariates = args.covariates
    phenotypes = args.phenotypes
    out_prefix = args.out_prefix
    response = args.response

    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    mt = io.import_table(input_path, input_type)

    if phenotypes:
        ht = hl.import_table(phenotypes,
                             types={'eid': hl.tstr, 'ukbb.centre' : hl.tstr},
                             missing=["", '""', "NA"],
                             impute=True,
                             force=True,
                             key='eid')
        mt = mt.annotate_cols(pheno=ht[mt.s])
    else:
        raise ValueError("param 'phenotypes' is not set!")


    # make one hot encoding for ukbb centre
    categories = mt.aggregate_cols(hl.agg.collect_as_set(mt.pheno['ukbb.centre']))
    categories = list(categories) 
    categories = [category for category in categories if category is not None]

    # Initialize a list to keep track of the new field names
    new_covariate_names = []
    
    # Annotate the MatrixTable with new fields for each category
    for category in categories:
        print(category)
        # Generate a safe variable name for the category
        encoded_var_name = f"ukbb.centre_{category}".replace(" ", "_")
        # Annotate the MatrixTable with a new column for this category
        mt = mt.annotate_cols(**{encoded_var_name: hl.if_else(mt.pheno['ukbb.centre'] == category, 1, 0)})
        
        # Keep track of this new variable name
        new_covariate_names.append(encoded_var_name)


    # Assuming 'split_covariates' contains the original list of covariate names including 'ukbb.centre'
    split_covariates = covariates.split(',')
    split_covariates = [cov for cov in split_covariates if cov != 'ukbb.centre']  # Remove 'ukbb.centre'
    split_covariates += new_covariate_names  # Add the new one-hot encoded variable names

    covariates = [mt.pheno[x] for x in split_covariates if x in mt.pheno]
    covariates += [mt[x] for x in split_covariates if x not in mt.pheno]  # Include one-hot encoded covariates
    covariates.insert(0, 1)  # Include intercept
    if response in list(mt.pheno):
        if min_maf_cutoff:
            mt = mt.filter_rows(variants.get_maf_expr(mt) > float(min_maf_cutoff))
        if mt.pheno[response].dtype == hl.dtype('float64'):
            reg = hl.linear_regression_rows(
                y=mt.pheno[response],
                x=mt.GT.n_alt_alleles(),
                covariates=covariates
            )
        else:
            raise TypeError("Response variable is not a float64!")
    else:
        raise ValueError("Response variable is not in phenotype file")

    # Get allele frequencies and annotate results
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='chromosomes to be parsed.')    
    parser.add_argument('--input_path', default=None, help='Path to input.')
    parser.add_argument('--input_type', default=None, help='Input type (vcf/mt/plink).')
    parser.add_argument('--min_maf_cutoff', default=None, help='Filter by minimum minor allele frequency')
    parser.add_argument('--response', default=None, help='Response variable for continuous traits')
    parser.add_argument('--covariates', default=None, help='List of covariates')
    parser.add_argument('--phenotypes', default=None, help='File path to phenotypes.')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)

