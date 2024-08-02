#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init

def DS_to_fake_GT(DS):
    """ Convert dosage to fake genotype """
    return (hl.case()
            .when(DS == 0, hl.parse_call("0/0"))
            .when(DS == 1, hl.parse_call("1/0"))
            .when(DS == 2, hl.parse_call("1/1"))
            .or_missing()
           )
    

def main(args):

    # parser
    input_path = args.input_path
    input_type = args.input_type
    covariates = args.covariates
    phenotypes = args.phenotypes
    pheno_file = args.pheno_file
    out_prefix = args.out_prefix

    hail_init.hail_bmrc_init_local('logs/hail/prefilter_phenotypes.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip

    main = io.import_table(input_path, input_type, calc_info = False)
    
    # setup phenotypes
    ht = hl.import_table(pheno_file,
                 types={'eid': hl.tstr},
                 missing=["",'""',"NA"],
                 impute=True,
                 force=True,
                 key='eid'
                 )

    # this is the table we keep track of
    main = main.annotate_cols(pheno=ht[main.s])
    main = main.annotate_rows(stats = hl.struct())

    # read in covariates
    phenotypes = phenotypes.strip().split(',')
    phenotypes = list(filter(None, phenotypes))
    split_covariates = covariates.strip().split(',')

    for phenotype in phenotypes:

        # point to original
        mt = main

        # need build covaraites in context of mt
        current_cols = [mt.pheno[x] for x in split_covariates]
        current_cols += [mt.pheno[phenotype]] 

        # filter columns based on phenotype availability
        mt = mt.filter_cols(
            hl.is_defined(current_cols)
        )

        # get allele count
        mt = mt.annotate_rows(
            AC = hl.agg.sum(mt.DS)
        )
        # annotate main matrix table with allele counts
        main = main.annotate_rows(
            stats = main.stats.annotate(
                AC = mt.rows()[main.row_key].AC
            )
        )
        main = main.annotate_rows(
            stats = main.stats.rename({'AC': "ac_" + str(phenotype)})
        )

    # export file
    main = main.rows().flatten()
    main.export(out_prefix + ".txt.gz")

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='Path to input.')
    parser.add_argument('--input_type', default=None, help='Input type (vcf/mt/plink).')
    parser.add_argument('--phenotypes', default=None, help='File path to phenotypes.')
    parser.add_argument('--pheno_file', default=None, help='File path to phenotypes.')
    parser.add_argument('--covariates', default=None, help='File path to phenotypes.')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)
