#!/usr/bin/env python3

import hail as hl
import argparse
import pandas
import os

from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import samples
from ko_utils import qc
from ko_utils import analysis

def main(args):
    
    # parser
    chrom = args.chrom
    input_phased_path = args.input_phased_path
    input_phased_type = args.input_phased_type
    input_unphased_path = args.input_unphased_path
    input_unphased_type = args.input_unphased_type
    input_imputed_path = args.input_imputed_path
    out_prefix = args.out_prefix
    out_type   = args.out_type
    vep_path   = args.vep_path
    final_variant_list = args.final_variant_list
    final_sample_list = args.final_sample_list    

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    
    # get tables
    mt1 = qc.get_table(input_path=input_phased_path, input_type=input_phased_type) # 12788
    mt2 = qc.get_table(input_path=input_unphased_path, input_type=input_unphased_type) # 11867 (for singletons)

    # remove withdrawn samples
    mt1 = samples.remove_withdrawn(mt1)
    mt2 = samples.remove_withdrawn(mt2)

    # filter by final samples
    if final_sample_list:
        ht_final_samples = hl.import_table(final_sample_list, no_header=True, key='f0', delimiter=',')
        mt1 = mt1.filter_cols(hl.is_defined(ht_final_samples[mt1.col_key]))
        mt2 = mt2.filter_cols(hl.is_defined(ht_final_samples[mt2.col_key]))
   
    # filter by final variants
    if final_variant_list:
        ht_final_variants = hl.import_table(final_variant_list, types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
        ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)
        mt1 = mt1.filter_rows(hl.is_defined(ht_final_variants[mt1.row_key]))
  
    # annotate with imputed data
    if input_imputed_path:
        imputed = hl.read_table(input_imputed_path)
        mt1 = mt1.annotate_cols(in_imputed = hl.is_defined(imputed[mt1.col_key].s)) 
        mt2 = mt2.annotate_cols(in_imputed = hl.is_defined(imputed[mt2.col_key].s))

    # Add QC fields that has been removed during phasing to mt1
    mt1 = mt1.annotate_entries(DP = mt2[(mt1.locus, mt1.alleles), mt1.s].DP)
    mt1 = mt1.annotate_entries(GQ = mt2[(mt1.locus, mt1.alleles), mt1.s].GQ)
 
    mt1 = mt1.annotate_cols(gq = hl.agg.stats(mt1.GQ), dp = hl.agg.stats(mt1.DP))
    mt1 = hl.sample_qc(mt1, name='sample_qc')
    mt1.cols().flatten().export(output=out_prefix + "_samples_phased.tsv.bgz")
 
    mt2 = mt2.annotate_cols(gq = hl.agg.stats(mt2.GQ), dp = hl.agg.stats(mt2.DP))
    mt2 = hl.sample_qc(mt2, name='sample_qc')
    mt2.cols().flatten().export(output=out_prefix + "_samples_unphased.tsv.bgz")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_phased_path', default=None, help='Path to input')
    parser.add_argument('--input_phased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_unphased_path', default=None, help='Path to input that contains singletons')
    parser.add_argument('--input_unphased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_imputed_path', default=None, help='Path to imputed data (GRCh38)')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--chrom', default=None, help='Chromosome to be used')
    parser.add_argument('--vep_path', default=None, help='path to a .vcf file containing annotated entries by locus and alleles')   
    parser.add_argument('--final_sample_list', default=None, help='Path to hail table with final samples to be included')
    parser.add_argument('--final_variant_list', default=None, help='Path to hail table with final variants to be included')     

    args = parser.parse_args()

    main(args)

