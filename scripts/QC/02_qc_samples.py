#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import samples
from ko_utils import qc

def main(args):
    
    # parser
    chrom = args.chrom
    input_phased_path = args.input_phased_path
    input_phased_type = args.input_phased_type
    input_unphased_path = args.input_unphased_path
    input_unphased_type = args.input_unphased_type
    out_prefix = args.out_prefix
    out_type   = args.out_type
    final_variant_list = args.final_variant_list
    final_sample_list = args.final_sample_list    
    combine_datasets = args.combine_datasets

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    
    # get tables
    if combine_datasets:
        mt1 = qc.get_table(input_path=input_phased_path, input_type=input_phased_type) 
        mt2 = qc.get_table(input_path=input_unphased_path, input_type=input_unphased_type) 
        mt = qc.union_phased_with_unphased(mt1, mt2)
    else:
        mt = qc.get_table(input_path=input_unphased_path, input_type=input_unphased_type)

    # remove withdrawn samples
    mt = samples.remove_withdrawn(mt)

    # filter by final samples
    if final_sample_list:
        ht_final_samples = hl.import_table(final_sample_list, no_header=True, key='f0', delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))
   
    # filter by final variants
    if final_variant_list:
        ht_final_variants = hl.import_table(final_variant_list, types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
        ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)
        mt = mt.filter_rows(hl.is_defined(ht_final_variants[mt.row_key]))
    
    # annotate with imputed data
    imp = genotypes.get_ukb_imputed_v3_bgen([chrom])
    mt = mt.annotate_cols(in_imputed = hl.is_defined(imp.cols()[mt.s]))
 
    mt = mt.annotate_cols(gq = hl.agg.stats(mt.GQ), dp = hl.agg.stats(mt.DP))
    mt = hl.sample_qc(mt, name='sample_qc')
    mt.cols().select('sample_qc','in_imputed').write(out_prefix + "_samples.ht", overwrite = True)
    mt.cols().flatten().export(out_prefix + "_samples.tsv.bgz")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_phased_path', default=None, help='Path to input')
    parser.add_argument('--input_phased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_unphased_path', default=None, help='Path to input that contains singletons')
    parser.add_argument('--input_unphased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--chrom', default=None, help='Chromosome to be used')
    parser.add_argument('--final_sample_list', default=None, help='Path to hail table with final samples to be included')
    parser.add_argument('--final_variant_list', default=None, help='Path to hail table with final variants to be included')     
    parser.add_argument('--combine_datasets', action='store_true', default=False, help='Combine phased data and unphased singletons into a single matrix table')     

    args = parser.parse_args()

    main(args)

