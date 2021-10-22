#!/usr/bin/env python3

import hail as hl
import argparse
from ukb_utils import hail_init
from ko_utils import qc

def main(args):
    
    # parser
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    final_variant_list = args.final_variant_list

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    mt = qc.get_table(input_path=input_path, input_type=input_type)
    
    if final_variant_list:
        ht_final_variants = hl.import_table(final_variant_list, types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
        ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)
        mt = mt.filter_rows(hl.is_defined(ht_final_variants[mt.row_key]))

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
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')


    args = parser.parse_args()

    main(args)

