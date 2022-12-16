#!/usr/bin/env python3


import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import variants


def main(args):

    out_prefix = args.out_prefix
    intervals = args.intervals
    extract = args.extract
    min_info = float(args.min_info)
    min_maf = float(args.min_maf)
    missing = float(args.missing)
    min_maf_by_case_control = args.min_maf_by_case_control
    pheno_file = args.pheno_file
    phenotype = args.phenotype
    trait = args.trait 
    checkpoint = args.checkpoint

    reference_genome = 'GRCh38'
    hail_init.hail_bmrc_init(log='logs/hail/filter_genotyoes.log', default_reference=reference_genome, min_block_size=128) 
    
    # read HailTable of intervals
    ht = hl.read_table(intervals)
    chromosomes = [x.replace('chr', '') for x in list(set(ht.contig.collect()))]
   
    # get imputation scores for chromosomes required. 
    mfis = list()
    for chrom in chromosomes:
        mfi = genotypes.get_ukb_imputed_v3_mfi(chrom)
        mfi = mfi.annotate(chrom = chrom)
        mfis.append(mfi)

    # get variants from MFI
    mfi = mfis[0].union(*mfis[1:])
    mfi = mfi.annotate(ref = hl.if_else(mfi.f6 == mfi.a1, mfi.a2, mfi.a1))
    mfi = mfi.annotate(variant = hl.delimit([hl.str(mfi.chrom), hl.str(mfi.position), mfi.ref, mfi.f6], ':'))
    mfi = mfi.key_by(**hl.parse_variant(mfi.variant,  reference_genome= 'GRCh37'))

    # get imputed genotypes
    mt = genotypes.get_ukb_imputed_v3_bgen(chromosomes)
    mt = mt.annotate_rows(info_score = mfi[mt.row_key].info)
    mt = mt.select_entries(mt.GT)

    # Filter to relevant samples
    if extract:
        ht_final_samples = hl.import_table(
            extract, no_header=True, key='f0', delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))

    # perform variant filtering after subsetting samples
    info_expr = mt.info_score > min_info
    maf_expr = variants.get_maf_expr(mt) > min_maf
    missing_expr = variants.get_missing_expr(mt) < missing
    mt = mt.filter_rows((info_expr) & (maf_expr) & (missing_expr))
    
    # perform liftover
    mt = variants.liftover(mt, fix_ref = False)
    
       # get intervals and filter file by the intervals
    hail_intervals = ht.intervals.collect()
    mt = hl.filter_intervals(mt, hail_intervals)
    
    # checkpoint after liftover
    if checkpoint:
        mt = mt.checkpoint(out_prefix + "_checkpoint.mt", overwrite = True)


    # add variant IDs
    mt = mt.annotate_rows(
                varid = hl.delimit(
                    [hl.str(mt.locus.contig),
                     hl.str(mt.locus.position),
                     mt.alleles[0],
                     mt.alleles[1]],
                    ':')
                ) 

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
        cases = mt.aggregate_cols(hl.agg.sum(mt[phenotype] == 1))
        controls = mt.aggregate_cols(hl.agg.sum(mt[phenotype] == 0)) 
        min_maf = hl.max(0.01, 25/(2 * hl.min([cases, controls]))).collect()[0]
        
        # write to outfile
        outfile = out_prefix + "_min_maf.tsv"
        with open(outfile, "w") as outfile:
            line = ("%s\t%d\t%d\t%f" % (phenotype, cases, controls, min_maf))  
            outfile.write(line + "\n")

    else:
        raise ValueError("param 'phenotypes' is not set! ")

    hl.export_vcf(mt, out_prefix + '.vcf.bgz')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--extract', default=None, help='Path to HailTable that contains the final samples included in the analysis.')
    parser.add_argument('--intervals', default=None, help='Path to HailTable that contains the significant genes from primary analysis.')
    parser.add_argument('--min_maf', default=0.01, help='What min_maf threshold should be used?')
    parser.add_argument('--min_info', default=0.8, help='What info threshold should be used?')
    parser.add_argument('--missing', default=0.1, help='What info threshold should be used?')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset (plink format)')
    parser.add_argument('--phenotype', default=None, help='Actual phenotpe for min_maf thresholding')
    parser.add_argument('--pheno_file', default=None, help='Path to phenotypes for min_maf thresholding')
    parser.add_argument('--min_maf_by_case_control', default=None, action="store_true", help='Should min_maf be set by case-control ratio?')
    parser.add_argument('--checkpoint', default=None, action="store_true", help='Should checkpoints be created along the way')
    parser.add_argument('--trait', default=None, help='What is the trait?')
    args = parser.parse_args()

    main(args)
