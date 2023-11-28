#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import samples
from ukb_utils import variants
from ko_utils import io

def main(args):
    
    in_file = args.in_file
    in_type = args.in_type
    input_annotation_path = args.input_annotation_path
    final_sample_list = args.final_sample_list
    final_variant_list = args.final_variant_list
    out_prefix = args.out_prefix
    out_type = args.out_type
    dbsnp_path = args.dbsnp_path
    annotate_snp_id = args.annotate_snp_id
    annotate_rsid = args.annotate_rsid

    hail_init.hail_bmrc_init_local('logs/hail/annotate_mt.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')

    mt = io.import_table(in_file, in_type, find_replace=(':-?nan', ':NaN'))
    mt = samples.remove_withdrawn(mt)

    if final_sample_list:
        ht_final_samples = hl.import_table(
                final_sample_list, 
                no_header=True, key='f0', 
                delimiter=',')
        mt = mt.filter_cols(
                hl.is_defined(ht_final_samples[mt.col_key]))

    if final_variant_list:
        ht_final_variants = hl.import_table(final_variant_list, 
                types={'locus': hl.tlocus(reference_genome='GRCh38'), 
                       'alleles': hl.tarray(hl.tstr)})
        ht_final_variants = ht_final_variants.key_by(
                ht_final_variants.locus, ht_final_variants.alleles)
        mt = mt.filter_rows(hl.is_defined(ht_final_variants[mt.row_key]))

    if input_annotation_path:
        consequence_annotations = hl.read_table(input_annotation_path)
        mt = mt.annotate_rows(consequence=consequence_annotations[mt.row_key]) 
    
    if annotate_rsid:
        if dbsnp_path:
            ht_dbsnp = variants.get_dbsnp_table(version=dbsnp_path, build='GRCh38')
            mt = mt.annotate_rows(rsid=ht_dbsnp.rows()[mt.row_key].rsid)
        else:
            raise ValueError("Missing argument 'dbsnp_path'")
        
    if annotate_snp_id:
        mt = mt.annotate_rows(
                varid = hl.delimit(
                    [hl.str(mt.locus.contig), 
                     hl.str(mt.locus.position), 
                     mt.alleles[0], 
                     mt.alleles[1]],
                    ':')
                )

    # always produce a matrix-table 
    if out_type not in "mt":
        io.export_table(mt, out_prefix, "mt")
    io.export_table(mt, out_prefix, out_type)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--in_file', default=None, required = True, help='Path to input')
    parser.add_argument('--in_type', default=None, required = True, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_annotation_path', default=None, help='Input for annotation path')
    parser.add_argument('--final_sample_list', default=None, help='Path to hail table with final samples to be included')
    parser.add_argument('--final_variant_list', default=None, help='Path to hail table with final variants to be included') 
    parser.add_argument('--dbsnp_path', default=None, help='what dbsnp version should be used for rsid annotation?') 
    parser.add_argument('--annotate_snp_id', default=None, action='store_true', help='Annotate variant id by chr:pos:a1:a2') 
    parser.add_argument('--annotate_rsid', default=None, action='store_true', help='Annotate variant id by chr:pos:a1:a2') 
    parser.add_argument('--out_prefix', default=None, required = True, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, required = True, help='Type of output dataset (options: mt, vcf, plink)')
 
    args = parser.parse_args()

    main(args)

