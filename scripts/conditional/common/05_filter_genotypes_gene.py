#!/usr/bin/env python3


import hail as hl
import pandas as pd
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

    print(f"running {phenotype} using Hail..")
    reference_genome = 'GRCh38'
    hail_init.hail_bmrc_init(log='logs/hail/filter_genotyoes.log', default_reference=reference_genome, min_block_size=128) 
    
    # get chromosomes and ranges to extract    
    ht = hl.read_table(intervals)
    df = ht.to_pandas()    
 
    outfile = out_prefix + "_intervals.txt"
    with open(outfile, "w") as outfile:
        for idx, row in df.iterrows():
            # iterate over rows
            gene = row['gene']
            contig = row['contig']
            start = row['start_with_padding']
            end = row['end_with_padding']
            interval = f"{contig}:{start}-{end}"
            hail_interval = hl.parse_locus_interval(interval, reference_genome='GRCh38')
            # create interval per region
            path = f"data/unphased/imputed/common_append_missing/ukb_imp_200k_common_append_missing_{contig}.mt"
            mt = hl.read_matrix_table(path)
            mt = hl.filter_intervals(mt, hl.literal([hail_interval])) 
            # create outfile 
            out_prefix_gene = out_prefix + "_" + gene + ".vcf.bgz"
            print(f"Exporting to {out_prefix_gene}.vcf.bgz")
            hl.export_vcf(mt, out_prefix_gene)
            line = ("%s\t%s\t%s\t%s\t%s\t%s" % (phenotype, gene, contig, start, end, out_prefix_gene))
            outfile.write(line + "\n")
    
    


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
