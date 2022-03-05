#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
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
    input_gnomad_path = args.input_gnomad_path
    input_imputed_path = args.input_imputed_path
    input_annotation_path = args.input_annotation_path
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
    
    # annotate with gnomAD
    if input_gnomad_path:
        gnomad_variants_ht = hl.import_vcf(input_gnomad_path, reference_genome ='GRCh38', force_bgz=True, array_elements_required=False).rows()
        mt = mt.annotate_rows(in_gnomad = hl.is_defined(gnomad_variants_ht[mt.row_key])) 
    
    # annotate with genebass
    if True:
        path = '/well/lindgren/UKBIOBANK/nbaya/resources/ukbb-exome-public/variant_results.mt'
        gb = hl.read_matrix_table(path)
        mt = mt.annotate_rows(in_genebass = hl.is_defined(gb.rows()[(mt.locus, mt.alleles)]))
    
    # annotate with imputed data
    if input_imputed_path:
        imputed = hl.read_table(input_imputed_path)
        mt = mt.annotate_rows(imputed_info = imputed[mt.row_key].info) 

    # add annotations from table
    consequence_annotations = hl.read_table(input_annotation_path)
    mt = mt.annotate_rows(consequence = consequence_annotations[mt.row_key]) 
   
    # get QC metrics
    mt = hl.variant_qc(mt, name='variant_qc')
    mt_rows = mt.rows().select('variant_qc','imputed_info', 'in_gnomad')
    mt_rows = mt_rows.annotate(worst_csq_for_variant_canonical = consequence_annotations[mt_rows.key].vep.worst_csq_for_variant_canonical)
    mt_rows.write(out_prefix + "_variants_qc.ht", overwrite=True)
    mt_rows.flatten().export(out_prefix + "_variants_qc.tsv.bgz")
 
    # annotate consequence category
    #category_annotations = analysis.annotation_case_builder(mt.consequence.vep.worst_csq_for_variant_canonical, use_loftee = True)
    #mt = mt.annotate_rows(consequence_category = category_annotations)
    
    # count URVs by samples 
    #mt_samples = analysis.count_urv_by_samples(mt)
    #mt_samples.entries().flatten().export(out_prefix + "_urv_by_samples.tsv.bgz")

    # count URVs by genes
    #mt = mt.explode_rows(mt.consequence.vep.worst_csq_by_gene_canonical)
    #mt_genes = analysis.count_urv_by_genes(mt, mt.consequence.vep.worst_csq_for_variant_canonical.gene_id)
    #mt_genes.entries().flatten().export(out_prefix + "_urv_by_genes.tsv.bgz")

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
    parser.add_argument('--final_sample_list', default=None, help='Path to hail table with final samples to be included')
    parser.add_argument('--final_variant_list', default=None, help='Path to hail table with final variants to be included')     
    parser.add_argument('--combine_datasets', action='store_true', default=False, help='Combine phased data and unphased singletons into a single matrix table')     


    args = parser.parse_args()

    main(args)

