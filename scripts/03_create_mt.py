#!/usr/bin/env python3

import hail as hl
import argparse

from gnomad.utils.vep import process_consequences
from ukb_utils import hail_init
from ukb_utils import samples
from ukb_utils import variants 
from ukb_utils import tables
from ko_utils import io
from ko_utils import variants as ko_utils_variants

def main(args):
    
    input_phased_path = args.input_phased_path
    input_phased_type = args.input_phased_type
    input_unphased_path = args.input_unphased_path
    input_unphased_type = args.input_unphased_type
    input_annotation_path = args.input_annotation_path
    final_sample_list = args.final_sample_list
    final_variant_list = args.final_variant_list
    out_prefix = args.out_prefix
    out_type = args.out_type
    varid = args.varid
    dbsnp = args.dbsnp

    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')

    phased = io.import_table(
            input_path=input_phased_path, 
            input_type=input_phased_type) # 12788
    
    unphased = io.import_table(
            input_path=input_unphased_path, 
            input_type=input_unphased_type) # 11867

    unphased = ko_utils_variants.filter_max_mac(unphased, 1)
    unphased = tables.order_cols(unphased, phased)
    mt = io.rbind_matrix_tables(unphased, phased)
    
    mt = samples.remove_withdrawn(mt)

    if final_sample_list:
        ht_final_samples = hl.import_table(
                final_sample_list, 
                no_header=True, key='f0', 
                delimiter=',')
        mt = mt.filter_cols(
                hl.is_defined(ht_final_samples[mt.col_key]))

    if final_variant_list:
        ht_final_variants = hl.import_table(
                final_variant_list, 
                types={'locus': hl.tlocus(reference_genome='GRCh38'), 
                       'alleles': hl.tarray(hl.tstr)})
        ht_final_variants = ht_final_variants.key_by(
                ht_final_variants.locus, ht_final_variants.alleles)
        mt = mt.filter_rows(hl.is_defined(ht_final_variants[mt.row_key]))

    if input_annotation_path:
        consequence_annotations = hl.read_table(input_annotation_path)
        mt = mt.annotate_rows(consequence=consequence_annotations[mt.row_key]) 
    
    if dbsnp:
        ht_dbsnp = variants.get_dbsnp_table(version=dbsnp, build='GRCh38')
        mt = mt.annotate_rows(rsid=ht_dbsnp.rows()[mt.row_key].rsid)
    
    if varid:
        mt = mt.annotate_rows(
                varid = hl.delimit(
                    [hl.str(mt.locus.contig), 
                     hl.str(mt.locus.position), 
                     mt.alleles[0], 
                     mt.alleles[1]],
                    ':'))

    io.export_table(mt, out_prefix, out_type)


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
    parser.add_argument('--dbsnp', default=None, help='what dbsnp version should be used for rsid annotation?') 
    parser.add_argument('--varid', default=None, action='store_true', help='Annotate variant id by chr:pos:a1:a2') 
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--chrom', default=None, help='Chromosome to be used')
    parser.add_argument('--vep_path', default=None, help='path to a .vcf file containing annotated entries by locus and alleles')   
 
    args = parser.parse_args()

    main(args)

