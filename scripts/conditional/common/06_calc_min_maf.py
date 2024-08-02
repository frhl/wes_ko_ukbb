#!/usr/bin/env python3


import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import variants


def main(args):

    out_prefix = args.out_prefix
    intervals = args.intervals
    min_maf_by_case_control = args.min_maf_by_case_control
    pheno_file = args.pheno_file
    phenotype = args.phenotype
    trait = args.trait 

    reference_genome = 'GRCh38'
    hail_init.hail_bmrc_init(log='logs/hail/filter_genotyoes.log', default_reference=reference_genome, min_block_size=128) 
    
    # get chromosomes and ranges to extract    
    ht = hl.read_table(intervals)
    chromosomes = [x.replace('chr', '') for x in list(set(ht.contig.collect()))]
    hail_intervals = ht.intervals.collect()
    
    # get all the relevant chromsomes
    mts = list()
    for chrom in chromosomes:
        path = f"data/unphased/imputed/common_append_missing/ukb_imp_200k_common_append_missing_chr{chrom}.mt"
        mt = hl.read_matrix_table(path)
        mts.append(mt) 
        print(f"loading {path}")
    
    # combine them by variants (across chroms)
    mt = mts[0].union_rows(*mts[1:])
    mt = hl.filter_intervals(mt, hail_intervals)
    
    # Import phenotypes for MAF thresholding
    if  pheno_file and min_maf_by_case_control and trait in "binary":
        ht = hl.import_table(pheno_file,
                     types={'eid': hl.tstr},
                     missing=["",'""',"NA"],
                     impute=True,
                     force=True,
                     key='eid')
        # subset min-maf by case controls
        mt = mt.annotate_cols(pheno=ht[mt.s][phenotype])
        cases = mt.aggregate_cols(hl.agg.sum(mt.pheno == True))
        controls = mt.aggregate_cols(hl.agg.sum(mt.pheno == False)) 
        min_maf = hl.max(0.01, 25/(2 * hl.min([cases, controls]))).collect()[0]
        
        # write to outfile
        outfile = out_prefix + "_min_maf.tsv"
        with open(outfile, "w") as outfile:
            line = ("%s\t%d\t%d\t%f" % (phenotype, cases, controls, min_maf))  
            outfile.write(line + "\n")

    else:
        raise ValueError("param 'phenotypes' is not set! ")

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--intervals', default=None, help='Path to HailTable that contains the significant genes from primary analysis.')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset (plink format)')
    parser.add_argument('--phenotype', default=None, help='Actual phenotpe for min_maf thresholding')
    parser.add_argument('--pheno_file', default=None, help='Path to phenotypes for min_maf thresholding')
    parser.add_argument('--min_maf_by_case_control', default=None, action="store_true", help='Should min_maf be set by case-control ratio?')
    parser.add_argument('--trait', default=None, help='What is the trait?')
    args = parser.parse_args()

    main(args)
