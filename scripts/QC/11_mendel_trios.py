#!/usr/bin/env python3

import hail as hl
import argparse
import pandas
import os

from gnomad.utils.vep import process_consequences
from ukb_utils import hail_init
from ukb_utils import samples
from ko_utils import qc
from ko_utils import analysis

def main(args):
    
    # parser
    #input_phased_path = args.input_phased_path
    #input_phased_type = args.input_phased_type
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
    mt = qc.get_table(input_path=input_unphased_path, input_type=input_unphased_type) # 11867 (for singletons

    # Remove withdrawn samples
    mt = samples.remove_withdrawn(mt)

    # filter by Duncan's final samples
    if final_sample_list:
        ht_final_samples = hl.import_table(final_sample_list, no_header=True, key='f0', delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))

    # filter by Ducan's final variants
    if final_variant_list:
        ht_final_variants = hl.import_table(final_variant_list, types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
        ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)
        mt = mt.filter_rows(hl.is_defined(ht_final_variants[mt.row_key]))

    # get mendel errors
    fam_path = samples.get_fam_path(app_id=11867,wes_200k_only=False,relateds_only=False)
    pedigree = hl.Pedigree.read(fam_path)
    all_errors, per_fam, per_sample, per_variant = hl.mendel_errors(mt['GT'], pedigree)
    mt = mt.annotate_rows(mendel=per_variant[mt.locus, mt.alleles])
    
    all_errors.export(out_prefix + '_mendel_all_errors.tsv.gz')
    per_variant.export(out_prefix + '_mendel_per_variant.tsv.gz')
    
    # Annotate with consequence
    consequence_annotations = hl.read_table(input_annotation_path)
    mt = mt.annotate_rows(consequence = consequence_annotations[mt.row_key]) 
    annotation = analysis.annotation_case_builder(mt.consequence.vep.worst_csq_for_variant_canonical, use_loftee = True)
    
    mt = mt.annotate_rows(csq_category = annotation) 
    mt = mt.annotate_rows(csq_variant = mt.consequence.vep.worst_csq_for_variant_canonical.most_severe_consequence)
    
    mt = mt.annotate_rows(
        count = hl.struct(
            coding_snp = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_snp(mt.alleles[0], mt.alleles[1])),
            coding_indel = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_indel(mt.alleles[0], mt.alleles[1]))
        )
    ) 
    
    # flatten and write
    mt_result = mt.select_rows('count','mendel', 'csq_category', 'csq_variant')
    mt_result.rows().flatten().export(out_prefix + "_mendel.tsv.gz")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    #parser.add_argument('--input_phased_path', default=None, help='Path to input')
    #parser.add_argument('--input_phased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
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


