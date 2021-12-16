#!/usr/bin/env python3

import hail as hl
import argparse
import pandas
import os

from gnomad.utils.vep import process_consequences
from ukb_utils import hail_init
from ukb_utils import samples
from ukb_utils import variants
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
    chrom   = args.chrom
    
    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') #
    
    # get tables
    mt1 = qc.get_table(input_path=input_phased_path, input_type=input_phased_type) # 12788
    mt2 = qc.get_table(input_path=input_unphased_path, input_type=input_unphased_type) # 11867 (for singletons

    # Remove withdrawn samples
    mt1 = samples.remove_withdrawn(mt1)
    mt2 = samples.remove_withdrawn(mt2)
    #mt2 = variants.liftover(mt2)

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

    # remove singletons since these can't be phased.
    mt1 = qc.filter_min_mac(mt1, 2)
    mt2 = qc.filter_min_mac(mt2, 2)
    
    # remove invariant sites
    mt1 = mt1.filter_rows(hl.agg.stats(mt1.GT.n_alt_alleles()).stdev == 0, keep=False)
    mt2 = mt2.filter_rows(hl.agg.stats(mt2.GT.n_alt_alleles()).stdev == 0, keep=False)

    # load pedigree
    fam_path = samples.get_fam_path(app_id=11867,wes_200k_only=False,relateds_only=False) # Files with True does not exist
    pedigree = hl.Pedigree.read(fam_path)
    
    # Remove invariant sites


    # setup trip matrix
    trio_dataset = hl.trio_matrix(mt2, pedigree, complete_trios=True)
    mt = hl.experimental.phase_trio_matrix_by_transmission(trio_dataset)
    mt = hl.experimental.explode_trio_matrix(mt)
   
    # clean up trio file
    mt = mt.drop('source_trio')
    mt = mt.select_entries('PBT_GT')
    mt = mt.annotate_entries(GT = mt.PBT_GT)
    mt = mt.select_entries('GT')
    mt = mt.drop('info')

    # Write VCF of trios with PBT_GT only
    sids = mt.s.collect()
    hl.export_vcf(mt, out_prefix + "_transmission.vcf")

    # shapeit4 phased data
    mt1 = mt1.filter_cols(hl.literal(set(sids)).contains(mt1.s))
    hl.export_vcf(mt1, out_prefix + "_estimation.vcf")

    # get call stats
    #info_field = 'info'
    #mt = mt.annotate_rows(**{info_field: hl.agg.call_stats(mt.PBT_GT, mt.alleles)})

    # annotate shapeit4 GTs with GTs phased by transmission and filter to hets
    #mt = mt.annotate_entries(GT_shapeit4 = mt1[mt.row_key, mt.col_key].GT)
    #mt = mt.filter_entries(hl.is_defined(mt.PBT_GT) & hl.is_defined(mt.GT_shapeit4))
    #mt = mt.filter_entries(mt.GT_shapeit4.is_het_ref())
    #mt = mt.filter_entries(mt.PBT_GT.is_het_ref())

    # annotate switch errors
    #mt = mt.annotate_rows(errors_sum=hl.agg.sum(mt.PBT_GT != mt.GT_shapeit4))
    #mt = mt.annotate_rows(hets_sum=hl.agg.sum(mt.PBT_GT.is_het_ref()))
    #mt = mt.annotate_rows(switch_errors=(mt.errors_sum / mt.hets_sum))

    # write to hail table
    #mt.rows().select('errors_sum','hets_sum','switch_errors','info').write(out_prefix + ".ht")
    
    #out = mt.select_entries(mt.PBT_GT, mt.GT, mt.GT_shapeit4)
    #out = out.drop(out.qual, out.filters)
    #out.entries().flatten().export(out_prefix + '_full.txt.bgz')

    
    # write to tab-seperated files
    #mt = mt.annotate_rows(AC1 = mt.info.AC[0])
    #mt = mt.annotate_rows(AC2 = mt.info.AC[1])
    #mt = mt.annotate_rows(AF1 = mt.info.AF[0])
    #mt = mt.annotate_rows(AF2 = mt.info.AF[1])
    #mt.rows().select('errors_sum','hets_sum','switch_errors','AC1','AC2','AF1','AF2').export(out_prefix + ".tsv.bgz")

    #mt.write(out_prefix + ".mt")

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


    # remove singletons since these can't be phased.
