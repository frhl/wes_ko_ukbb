#!/usr/bin/env python3

import hail as hl
import argparse
import pandas
import os

from gnomad.utils.vep import process_consequences
from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import samples
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

    # Remove withdrawn samples
    mt1 = samples.remove_withdrawn(mt1)
    mt2 = samples.remove_withdrawn(mt2)

    # filter by Duncan's final samples
    if final_sample_list:
        ht_final_samples = hl.import_table(final_sample_list, no_header=True, key='f0', delimiter=',')
        mt1 = mt1.filter_cols(hl.is_defined(ht_final_samples[mt1.col_key]))
        mt2 = mt2.filter_cols(hl.is_defined(ht_final_samples[mt2.col_key]))
   
    # filter by Ducan's final variants
    if final_variant_list:
        ht_final_variants = hl.import_table(final_variant_list, types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
        ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)
        mt1 = mt1.filter_rows(hl.is_defined(ht_final_variants[mt1.row_key]))
        mt2 = mt2.filter_rows(hl.is_defined(ht_final_variants[mt2.row_key]))

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

    # counts
    n_mt1 = mt1.count()
    n_mt2 = mt2.count()
    print(f'chr{chrom}: phased data count: {n_mt1}')
    print(f'chr{chrom}: unphased (singleton) data count: {n_mt2}')


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

