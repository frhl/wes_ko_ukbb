#!/usr/bin/env python3


import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import variants
from ukb_utils import samples
from ko_utils import io
from ko_utils.variants import filter_min_maf, filter_missing


def main(args):

    out_prefix = args.out_prefix
    gene_table = args.gene_table
    extract = args.extract
    min_info = float(args.min_info)
    min_maf = float(args.min_maf)
    missing = float(args.missing)
    padding = int(float(args.padding))
    min_maf_by_case_control = args.min_maf_by_case_control
    pheno_file = args.pheno_file
    phenotype = args.phenotype
    trait = args.trait 

    reference_genome = 'GRCh38'
    hail_init.hail_bmrc_init_local(
        'logs/hail/hail_conditional_tables.log',
        reference_genome)
    hl._set_flags(no_whole_stage_codegen='1')

    # import table of start/end contigs
    ht = hl.import_table(
        gene_table,
        force=True,
        delimiter='\t',
        missing='',
        types={
            'contig': hl.tstr,
            'start': hl.tint,
            'end': hl.tint})
    ht = ht.annotate(contig=hl.delimit(['chr', ht.contig], ''))
    
    # get refernece genome with lengths
    rg = hl.get_reference('GRCh38')
    chr_lens = hl.literal(rg.lengths)
    ht = ht.annotate(contig_length=chr_lens[ht.contig])

    # annotate regions with padding
    ht = ht.annotate(ranges=hl.array([ht.start, ht.end]))
    ht = ht.annotate(
        end=hl.max(ht.ranges),
        start=hl.min(ht.ranges))
    ht = ht.annotate(
        end_with_padding=hl.min(hl.array([ht.end + padding, ht.contig_length])),
        start_with_padding=hl.max(hl.array([ht.start - padding, 1])))
    ht = ht.drop(ht.ranges)

    # check if locus is valid, i.e. does it go above contig
    # length? This is only an issue when padding is large
    ht = ht.annotate(
        start_valid=hl.is_valid_locus(
            ht.contig,
            ht.start_with_padding,
            reference_genome))
    ht = ht.annotate(
        end_valid=hl.is_valid_locus(
            ht.contig,
            ht.end_with_padding,
            reference_genome))
    ht = ht.annotate(
            valid_intervals=ht.start_valid & ht.end_valid)

    defined_coords = hl.is_defined(ht.valid_intervals)
    n_drop = ht.filter(~defined_coords).count()
    if n_drop > 0:
        print(f'Dropping {n_drop} undefined start/end genomic coordinates.')
    ht = ht.filter(defined_coords)
    invalid_intervals = ht.filter(~ht.valid_intervals).count()
    assert invalid_intervals == 0, 'Some of the supplied intervals are outside current chromosome contig'
    
    # annotate intervals
    ht = ht.annotate(
        intervals=hl.locus_interval(
            ht.contig,
            ht.start_with_padding,
            ht.end_with_padding,
            True,
            True,
            reference_genome=reference_genome))

    
    chromosomes = [x.replace('chr', '') for x in list(set(ht.contig.collect()))]
    
    # Get imputation scores
    mfis = list()
    for chrom in chromosomes:
        mfi = genotypes.get_ukb_imputed_v3_mfi(chrom)
        mfi = mfi.annotate(chrom = chrom)
        mfis.append(mfi)
    
    #hts = [genotypes.get_ukb_imputed_v3_mfi(chrom) for chrom in chromosomes]
    mfi = mfis[0].union(*mfis[1:])
    mfi = mfi.annotate(ref = hl.if_else(mfi.f6 == mfi.a1, mfi.a2, mfi.a1))
    mfi = mfi.annotate(variant = hl.delimit([hl.str(mfi.chrom), hl.str(mfi.position), mfi.ref, mfi.f6], ':'))
    mfi = mfi.key_by(**hl.parse_variant(mfi.variant,  reference_genome= 'GRCh37'))
    mfi = mfi.filter(mfi.info >= min_info)
    
    # get imputed genotypes
    mt = genotypes.get_ukb_imputed_v3_bgen(chromosomes)
    mt = samples.remove_withdrawn(mt)
    mt = mt.filter_rows(hl.is_defined(mfi[mt.row_key]))
    mt = mt.select_entries(mt.GT)

    # Filter samples
    if extract:
        ht_final_samples = hl.import_table(
            extract, no_header=True, key='f0', delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))

    # filter variants based on missingness, 
    # drop info field so that we don't have to recalculate
    mt = io.recalc_info(mt)
    mt = filter_missing(mt, missing)
    mt = filter_min_maf(mt, min_maf)
    mt = mt.drop('info')

    # perform liftover to GRCh38 and filter to intervals
    mt = variants.liftover(mt)
    hail_intervals = ht.intervals.collect()
    mt = hl.filter_intervals(mt, hail_intervals)

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
        mt = mt.annotate_cols(pheno=ht[mt.s])
        cases = mt.aggregate_cols(hl.agg.sum(mt.pheno[phenotype] == 1))
        controls = mt.aggregate_cols(hl.agg.sum(mt.pheno[phenotype] == 0)) 
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
    parser.add_argument('--padding', default=0, help='How much extra padding should be included around genes (upstream/downsream)?')
    parser.add_argument('--extract', default=None, help='Path to HailTable that contains the final samples included in the analysis.')
    parser.add_argument('--gene_table', default=None, help='Path to HailTable that contains the significant genes from primary analysis.')
    parser.add_argument('--min_maf', default=0.01, help='What min_maf threshold should be used?')
    parser.add_argument('--min_info', default=0.8, help='What info threshold should be used?')
    parser.add_argument('--missing', default=0.1, help='What info threshold should be used?')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset (plink format)')
    parser.add_argument('--phenotype', default=None, help='Actual phenotpe for min_maf thresholding')
    parser.add_argument('--pheno_file', default=None, help='Path to phenotypes for min_maf thresholding')
    parser.add_argument('--min_maf_by_case_control', default=None, action="store_true", help='Should min_maf be set by case-control ratio?')
    parser.add_argument('--trait', default=None, help='What is the trait?')
    args = parser.parse_args()

    main(args)
