#!/usr/bin/env python3

import hail as hl
import argparse
import pandas
import os

from gnomad.utils.vep import process_consequences
from ukb_utils import hail_init
from ukb_utils import genotypes
from ko_utils import qc
from ko_utils import analysis


def summary_count_urv(mt):
    annotation_by_variant = analysis.annotation_case_builder(mt.vep.worst_csq_for_variant_canonical, mt.dbnsfp)
    mt = mt.annotate_rows(consequence_category = annotation_by_variant)
    mt = analysis.counts_urv_by_samples(mt)
    return mt.cols()

def summary_count_homozygous_urv(mt):
    annotation_by_variant = analysis.annotation_case_builder(mt.vep.worst_csq_for_variant_canonical, mt.dbnsfp)
    mt = mt.annotate_rows(consequence_category = annotation_by_variant)
    mt = analysis.count_homozygous_urv_by_samples(mt)
    return mt.cols()


def main(args):
    
    # parser
    chrom = args.chrom
    input_phased_path = args.input_phased_path
    input_phased_type = args.input_phased_type
    input_unphased_path = args.input_unphased_path
    input_unphased_type = args.input_unphased_type
    input_gnomad_path = args.input_gnomad_path
    input_imputed_path = args.input_imputed_path
    input_annotation_path = args.input_annotation_path
    out_prefix = args.out_prefix
    out_type   = args.out_type
    vep_path   = args.vep_path
    final_variant_list = args.final_variant_list
    final_sample_list = args.final_sample_list    

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    #hl._set_flags(no_whole_stage_codegen='1') # from zulip
    
    # get tables
    mt1 = qc.get_table(input_path=input_phased_path, input_type=input_phased_type) # 12788
    mt2 = qc.get_table(input_path=input_unphased_path, input_type=input_unphased_type) # 11867 (for singletons)

    # filter by final samples
    if final_sample_list:
        ht_final_samples = hl.import_table(final_sample_list, no_header=True, key='f0')
        mt1 = mt1.filter_cols(hl.is_defined(ht_final_samples[mt1.col_key]))
        mt2 = mt2.filter_cols(hl.is_defined(ht_final_samples[mt2.col_key]))
        print(f'chr{chrom}: using final samples..')
   
    # filter by final variants
    if final_variant_list:
        ht_final_variants = hl.import_table(final_variant_list, types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
        ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)
        mt1 = mt1.filter_rows(hl.is_defined(ht_final_variants[mt1.row_key]))
        mt2 = mt2.filter_rows(hl.is_defined(ht_final_variants[mt2.row_key]))
        print(f'chr{chrom}: using final variants..')
    
    # annotate with gnomAD
    if input_gnomad_path:
        gnomad_variants_ht = hl.import_vcf(input_gnomad_path, reference_genome ='GRCh38', force_bgz=True, array_elements_required=False).rows()
        mt1 = mt1.annotate_rows(inGnomAD = hl.is_defined(gnomad_variants_ht[mt1.row_key]))
        mt2 = mt2.annotate_rows(inGnomAD = hl.is_defined(gnomad_variants_ht[mt2.row_key]))
        print(f'chr{chrom}: annotating with gnomAD..')

    # annotate with imputed data
    if input_imputed_path:
        imputed = hl.read_table(input_imputed_path)
        mt1 = mt1.annotate_rows(imputed_info = imputed[mt1.row_key].info) 
        mt2 = mt2.annotate_rows(imputed_info = imputed[mt2.row_key].info)
        print(f'chr{chrom}: annotating with imputed INFO score..')
    
    # Add QC fields that has been removed during phasing to mt1
    mt1 = mt1.annotate_entries(DP = mt2[(mt1.locus, mt1.alleles), mt1.s].DP)
    mt1 = mt1.annotate_entries(GQ = mt2[(mt1.locus, mt1.alleles), mt1.s].GQ)
 
    # save sample stats
    print(f'chr{chrom}: Writing out sample stats to {out_prefix}_samples*') 
    mt1 = mt1.annotate_cols(gq = hl.agg.stats(mt1.GQ), dp = hl.agg.stats(mt1.DP))
    mt1 = hl.sample_qc(mt1, name='sample_qc')
    mt1.cols().select('sample_qc', 'gq', 'dp').flatten().export(output=out_prefix + "_samples_phased.tsv.bgz")
 
    mt2 = mt2.annotate_cols(gq = hl.agg.stats(mt2.GQ), dp = hl.agg.stats(mt2.DP))
    mt2 = hl.sample_qc(mt2, name='sample_qc')
    mt2.cols().select('sample_qc', 'gq', 'dp').flatten().export(output=out_prefix + "_samples_unphased.tsv.bgz")
    
    # add annotations from table
    consequence_annotations = hl.read_table(input_annotation_path)
    mt1 = mt1.annotate_rows(consequence = consequence_annotations[mt1.row_key]) 
    mt2 = mt2.annotate_rows(consequence = consequence_annotations[mt2.row_key]) 
    
    # write out variant stats
    print(f'chr{chrom}: Writing out variants stats to {out_prefix}_variants*')
    mt1 = hl.variant_qc(mt1, name='variant_qc')
    ht1_rows_filter = mt1.rows()
    ht1_rows_filter.write(out_prefix + "_variants_phased.ht", overwrite=True)
 
    mt2 = hl.variant_qc(mt2, name='variant_qc')
    ht2_rows_filter = mt2.rows()
    ht2_rows_filter.write(out_prefix + "_variants_unphased.ht", overwrite=True)

    # write out summary stats
    print(f"chr{chrom}: writing out variant summaries")
    summary_count_urv(mt1).export(out_prefix + "variants_summary_phased.tsv.bgz")
    summary_count_urv(mt2).export(out_prefix + "variants_summary_unphased.tsv.bgz")
    
    # get homozygous stats
    summary_count_homozygous_urv(mt1).export(out_prefix + "variants_homozygous_summary_phased.tsv.bgz")
    summary_count_homozygous_urv(mt2).export(out_prefix + "variants_homozygous_summary_unphased.tsv.bgz")


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_phased_path', default=None, help='Path to input')
    parser.add_argument('--input_phased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_unphased_path', default=None, help='Path to input that contains singletons')
    parser.add_argument('--input_unphased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_gnomad_path', default=None, help='Get path to gnomAD (GRCh38)')
    parser.add_argument('--input_imputed_path', default=None, help='Path to imputed data (GRCh38)')
    parser.add_argument('--input_annotation_path', default=None, help='path to HailTable with VEP and dbNSFP annotations')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--chrom', default=None, help='Chromosome to be used')
    parser.add_argument('--vep_path', default=None, help='path to a .vcf file containing annotated entries by locus and alleles')   
    parser.add_argument('--final_sample_list', default=None, help='Path to hail table with final samples to be included')
    parser.add_argument('--final_variant_list', default=None, help='Path to hail table with final variants to be included')     

    args = parser.parse_args()

    main(args)

