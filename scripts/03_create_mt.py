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

def main(args):
    
    # parser
    input_phased_path = args.input_phased_path
    input_phased_type = args.input_phased_type
    input_unphased_path = args.input_unphased_path
    input_unphased_type = args.input_unphased_type
    input_annotation_path = args.input_annotation_path
    final_sample_list = args.final_sample_list
    final_variant_list = args.final_variant_list
    out_prefix = args.out_prefix
    out_type   = args.out_type
    vep_path   = args.vep_path
    
    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') #
    
    # get tables
    mt1 = qc.get_table(input_path=input_phased_path, input_type=input_phased_type) # 12788
    mt2 = qc.get_table(input_path=input_unphased_path, input_type=input_unphased_type) # 11867 (for singletons

    # filter by Duncan's final samples
    if final_sample_list:
        ht_final_samples = hl.import_table(FINAL_SAMPLE_LIST, no_header=True, key='f0')
        mt1 = mt1.filter_cols(hl.is_defined(ht_final_samples[mt1.col_key]))
        mt2 = mt2.filter_cols(hl.is_defined(ht_final_samples[mt2.col_key]))
   
    # filter by Ducan's final variants
    if final_variant_list:
        ht_final_variants = hl.import_table(FINAL_VARIANT_LIST, types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
        ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)
        mt1 = mt1.filter_rows(hl.is_defined(ht_final_variants[mt1.row_key]))
        mt2 = mt2.filter_rows(hl.is_defined(ht_final_variants[mt2.row_key]))

    # Add QC fields that has been removed during phasing to mt1
    mt1 = mt1.annotate_entries(DP = mt2[(mt1.locus, mt1.alleles), mt1.s].DP)
    mt1 = mt1.annotate_entries(GQ = mt2[(mt1.locus, mt1.alleles), mt1.s].GQ)

    # note: need to translate ids to combine later!
    mt1 = qc.annotate_european(mt1)
    mt2 = qc.annotate_european(mt2)

    ### Variant filtering/annotations
    # Using mt2 as a singleton matrix so remove those with AC > 1
    mt2 = qc.filter_max_mac(mt2, 1)
    mt2 = qc.filter_min_mac(mt2, 1)

    # add annotations from table
    consequence_annotations = hl.read_table(input_annotation_path)
    mt1 = mt1.annotate_rows(consequence = consequence_annotations[mt1.row_key]) 
    mt2 = mt2.annotate_rows(consequence = consequence_annotations[mt2.row_key]) 

    # Run VEP + gnomAD variant annotations
    #mt1 = process_consequences(hl.vep(mt1, "utils/configs/vep_env.json"))
    #mt2 = process_consequences(hl.vep(mt2, "utils/configs/vep_env.json"))
    
    # Annotate with REVEL+CADD scores
    #mt1 = analysis.annotate_dbnsfp(mt1, vep_path)
    #mt2 = analysis.annotate_dbnsfp(mt2, vep_path)
   
    # Annotate for each consequence
    #mt1 = mt1.explode_rows(mt1.vep.worst_csq_by_gene_canonical)
    #mt2 = mt2.explode_rows(mt2.vep.worst_csq_by_gene_canonical)

    # annotate consequnece categories 
    #mt1 = analysis.variant_csqs_category_builder(mt1)
    #mt2 = analysis.variant_csqs_category_builder(mt2)
        
    # annotate Gene (In the future just use vep.worst_csq_by_gene_canonical downstream..) 
    #mt1 = mt1.annotate_rows(vep = mt1.vep.annotate(Gene = mt1.vep.worst_csq_by_gene_canonical.gene_id))
    #mt2 = mt2.annotate_rows(vep = mt2.vep.annotate(Gene = mt2.vep.worst_csq_by_gene_canonical.gene_id))
    
    # By default add snpid id annotation
    mt1 = qc.annotate_snpid(mt1)
    mt2 = qc.annotate_snpid(mt2)
    
    # annotate rsid
    mt1 = qc.annotate_rsid(mt1)
    mt2 = qc.annotate_rsid(mt2)
    mt1 = qc.default_to_snpid_when_missing_rsid(mt1)
    mt2 = qc.default_to_snpid_when_missing_rsid(mt2)
    
    # export files
    mt1.write(out_prefix + ".mt")
    mt2.write(out_prefix + "_singletons.mt")


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_phased_path', default=None, help='Path to input')
    parser.add_argument('--input_phased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_unphased_path', default=None, help='Path to input that contains singletons')
    parser.add_argument('--input_unphased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_annotation_path', default=None, help='Input for annotation path')
    parser.add_argument('--final_sample_list', default=None, help='Path to hail table with final samples to be included')
    parser.add_argument('--final_variant_list', default=None, help='Path to hail table with final variants to be included') 
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--chrom', default=None, help='Chromosome to be used')
    parser.add_argument('--vep_path', default=None, help='path to a .vcf file containing annotated entries by locus and alleles')   
 
    args = parser.parse_args()

    main(args)

