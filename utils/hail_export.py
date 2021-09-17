#!/usr/bin/env python3

#' @description Compound heterozygos HAIL pipeline
#' @todo Integrate SAIGE-GENE+
#' @DONE @todo Convert MatrixTable (sample, variant)-pairs to 'long' format
#' @DONE @todo Long format table should also contain variants found and their strand.
#' @DONE @todo Test burden_dosage function
#' @DONE @todo Create function that subsets variants of MODERATE impact
#' @DONE @todo hardy eq test for each variant
#' @DONE @todo add mt.repartition to different commands
#' @DONE @todo check worst consequence
#' @DONE @todo export ultra rare variants
#' @todo check which allele is included in VEP? What about MAF thresholds?


import hail as hl
import argparse
import pandas
import os

from ukb_utils import hail_init
from ukb_utils import genotypes
from ko_utils import qc
from ko_utils import analysis

def main(args):
    
    # parser
    input_phased_path = args.input_phased_path
    input_phased_type = args.input_phased_type
    input_unphased_path = args.input_unphased_path
    input_unphased_type = args.input_unphased_type
    out_prefix = args.out_prefix
    out_type   = args.out_type
    vep_path   = args.vep_path
    
    chrom      = int(args.chrom)
    maf_max    = (args.maf_max)
    maf_min    = (args.maf_min)
    hwe        = (args.hwe)
    missing    = args.missing
    
    get_related = args.get_related
    get_unrelated = args.get_unrelated
    get_europeans = args.get_europeans
    vep_filter = args.vep_filter
    export_ko_dosage_matrix = args.export_ko_dosage_matrix
    #export_ko_samples = args.export_ko_samples
    export_burden = args.export_burden
    export_ko_probability = args.export_ko_probability
    export_fake_vcf = args.export_fake_vcf

    # run parser
    hail_init.hail_bmrc_init('logs/hail/hail_export.log', 'GRCh38')
    mt1 = qc.get_table(input_path=input_phased_path, input_type=input_phased_type) # 12788
    mt2 = qc.get_table(input_path=input_unphased_path, input_type=input_unphased_type) # 11867 (for singletons)

    ### Sample filtering
    if get_related and not get_unrelated:
        mt1 = qc.filter_to_unrelated(mt1, get_related = True)
        mt2 = qc.filter_to_unrelated(mt2, get_related = True)

    if get_unrelated and not get_related:
        mt1 = qc.filter_to_unrelated(mt1, get_related = False)
        mt2 = qc.filter_to_unrelated(mt2, get_related = False)

	# note: need to translate ids to combine later!
    mt1 = qc.translate_sample_ids(mt1, 12788, 11867)
    if get_europeans:
        mt1 = qc.filter_to_european(mt1)
        mt2 = qc.filter_to_european(mt2)

    ### Variant filtering/annotations
    # Using mt2 as a singleton refereence, so remove those with AC > 1
    mt2 = qc.filter_max_mac(mt2, 1)

    if missing:
        mt1 = qc.filter_min_missing(mt1, float(missing))
        mt2 = qc.filter_min_missing(mt2, float(missing))

    if maf_max:
        mt1 = qc.filter_max_maf(mt1, float(maf_max))

    if maf_min:
        mt1 = qc.filter_min_maf(mt1, float(maf_min))

    if hwe:
        mt1 = qc.filter_hwe(mt1, float(hwe))

    if vep_path:
        mt1 = analysis.annotate_vep(mt1, vep_path)
        mt2 = analysis.annotate_vep(mt2, vep_path)
        if vep_filter: 
            mt1 = analysis.filter_vep(mt1, 'consequence_category', vep_filter)
            mt2 = analysis.filter_vep(mt2, 'consequence_category', vep_filter) 

    #### get stats

    if export_burden:

        # Count burden per gene per individual
        mt1_cat = analysis.gene_burden_category_annotations_per_sample(mt1)
        mt2_cat = analysis.gene_burden_category_annotations_per_sample(mt2)

        # combine singleton table and full table
        res = mt1_cat.annotate_entries(singletons = mt2_cat[(mt1_cat.Gene, mt1_cat.consequence_category), mt1_cat.s].n)
        res = res.annotate_entries(singletons = hl.if_else(hl.is_missing(res.singletons),0,res.singletons))
        res = res.annotate_entries(total = res.n + res.singletons)
        res = res.entries()
        res = res.filter(res.total > 0)

        # export data
        res.export(out_prefix + '_burden.tsv.gz')

    if export_ko_probability:

        # determine probability of being a ko
        mt_ko = analysis.gene_csqs_calc_pKO(mt1, mt2, 'dosage')
        mt_ko_entries = mt_ko.entries()
        mt_ko_entries = mt_ko_entries.filter(mt_ko_entries.pKO>0)

        # export data
        mt_ko_entries.export(out_prefix + '_ko_prob.tsv.gz')

    if export_ko_dosage_matrix:
        # ignores singletons (but gives an overview)
        mt_ko_matrix = analysis.gene_csqs_case_builder(mt1)
        mt_ko_matrix.export(out_prefix + '_ko_matrix_no_singletons.tsv.gz')

    if export_fake_vcf:
        out = analysis.gene_csqs_calc_pKO_pseudoSNP(mt1, mt2, chrom)
        qc.export_table(out, out_prefix = out_prefix + "_ko", out_type = 'vcf')
        qc.export_table(out, out_prefix = out_prefix + "_ko", out_type = 'plink')

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_phased_path', default=None, help='Path to input')
    parser.add_argument('--input_phased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_unphased_path', default=None, help='Path to input that contains singletons')
    parser.add_argument('--input_unphased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    # filtering variants
    parser.add_argument('--chrom', default=None, help='Chromosome to be used')
    parser.add_argument('--maf_max', default=None, help='Select all variants with a maf less than the indicated value')
    parser.add_argument('--maf_min', default=None, help='Select all variants with a maf greater than the indicated values')
    parser.add_argument('--hwe', default=None, help='Filter variants by HWE threshold')
    parser.add_argument('--missing', default=0.05, help='Filter variants by missingness threshold')
    # filtering samples
    parser.add_argument('--get_related', action='store_true', help='Select all samples that are related')
    parser.add_argument('--get_unrelated', action='store_true', help='Select all samples that are unrelated')
    parser.add_argument('--get_europeans', action='store_true', help='Filter to genetically confimed europeans?')
    # out
    parser.add_argument('--export_ko_probability', action='store_true', help='Exports the KO probability.')
    parser.add_argument('--export_burden', action='store_true', help='Export burden variant count by gene and and individuals.')
    parser.add_argument('--export_fake_vcf', action='store_true', help='Export a "fake" VCF file that contains KO probabilities as DP field..')
    parser.add_argument('--vep_path', default=None, help='path to a .vcf file containing annotated entries by locus and alleles')
    parser.add_argument('--vep_filter', nargs='+', help='Filter consequence_category by mutations e.g., "damaging_missense" or "ptv"')
    parser.add_argument('--export_ko_dosage_matrix', action='store_true', help='Generate a gene x sample matrix with KO status')
    
    args = parser.parse_args()

    main(args)

