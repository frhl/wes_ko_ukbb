#!/usr/bin/env python3

import hail as hl
import argparse
from ukb_utils import hail_init
from ukb_utils import samples
from ko_utils import qc

def main(args):
    
    # parser
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    final_variant_list = args.final_variant_list
    final_sample_list = args.final_sample_list
    input_gnomad_path = args.input_gnomad_path
    input_imputed_path = args.input_imputed_path
    
    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/annotate_genebass_gnomad_imputed.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    mt = qc.get_table(input_path=input_path, input_type=input_type)
    mt = samples.remove_withdrawn(mt)

    if final_variant_list:
        ht_final_variants = hl.import_table(final_variant_list, types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
        ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)
        mt = mt.filter_rows(hl.is_defined(ht_final_variants[mt.row_key]))

    if final_sample_list:
        ht_final_samples = hl.import_table(final_sample_list, no_header=True, key='f0', delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))
   
    if input_gnomad_path:
        gnomad_variants_ht = hl.import_vcf(input_gnomad_path, reference_genome ='GRCh38', force_bgz=True, array_elements_required=False).rows()
        mt = mt.annotate_rows(inGnomAD = hl.is_defined(gnomad_variants_ht[mt.row_key]))

    if input_imputed_path:
        imputed = hl.read_table(input_imputed_path)
        mt = mt.annotate_rows(imputed_info = imputed[mt.row_key].info)  
    
    # do variant QC
    mt = hl.variant_qc(mt, name = 'variant_qc')

    # get variant level data 
    path = '/well/lindgren/UKBIOBANK/nbaya/resources/ukbb-exome-public/variant_results.mt'
    gb = hl.read_matrix_table(path)
    mt = mt.annotate_rows(gb_AC = gb.rows()[(mt.locus, mt.alleles)].AC)
    mt = mt.annotate_rows(gb_AF = gb.rows()[(mt.locus, mt.alleles)].AF)
    mt = mt.annotate_rows(gb_gene = gb.rows()[(mt.locus, mt.alleles)].gene)
    mt = mt.annotate_rows(gb_annotation = gb.rows()[(mt.locus, mt.alleles)].annotation)
    
    # write file
    ht = mt.rows()
    ht.flatten().export(out_prefix + '.tsv.bgz')

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--final_variant_list', default=None, help='Path to hail table with final variants to be included')
    parser.add_argument('--final_sample_list', default=None, help='Path to hail table with final variants to be included')
    parser.add_argument('--input_gnomad_path', default=None, help='Path to input')
    parser.add_argument('--input_imputed_path', default=None, help='Path to imputed data')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')


    args = parser.parse_args()

    main(args)

