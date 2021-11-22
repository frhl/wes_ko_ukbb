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

    # get trios only
    ht = samples.get_fam(app_id=11867, wes_200k_only=False)
    trios = ht.filter((ht.PAT != '0') & (ht.MAT != '0'))

    # Filter the two matrix tables to the trios defined
    #mt1 = mt1.filter_cols(hl.is_defined(trios[mt1.s].FID))
    #mt2 = mt2.filter_cols(hl.is_defined(trios[mt2.s].FID))
    #imt1_n_trios = mt1.count()[1]
    #mt2_n_trios = mt2.count()[1]
    #assert mt1_n_trios == mt2_n_trios
    #print(f"chr{chrom}: Filtering to {mt1_n_trios} trios")

    # load pedigree
    fam_path = samples.get_fam_path(app_id=11867,wes_200k_only=False,relateds_only=False) # Files with True does not exist
    pedigree = hl.Pedigree.read(fam_path)
    print(f'chr{chrom}: loaded pedigree')

    # phase by transmission
    trio_dataset_mt2 = hl.trio_matrix(mt2, pedigree, complete_trios=True)
    n = trio_dataset_mt2.count()
    print(f"trio count {n} (mt2)")
    
    trio_dataset_mt1 = hl.trio_matrix(mt1, pedigree, complete_trios=True)
    n = trio_dataset_mt1.count()
    print(f"trio count {n} (mt1)")
    
    mt = hl.experimental.phase_trio_matrix_by_transmission(trio_dataset_mt2)
    n = mt.count()
    print(f'chr{chrom}: phased trios with count of {n}')

    # annotate shapeit4 GTs with GTs phased by transmission
    mt = mt.annotate_entries(GT_shapeit4 = mt1[mt.row_key, mt.col_key].GT)
    
    mt = mt.annotate_rows(switch_errors_count=hl.agg.sum(mt.proband_entry.PBT_GT != mt.GT_shapeit4))
    mt = mt.annotate_rows(switch_errors=hl.agg.fraction(mt.proband_entry.PBT_GT != mt.GT_shapeit4))
    #print('chr{chrom}: annotated rows')

    # write to file
    mt.rows().select('switch_errors','switch_errors_count').export(out_prefix + ".tsv.bgz")
    SE = mt.filter_rows(mt.switch_errors > 0).count()
    print(f"chr{chrom}: Rows with switch errors count: {SE}")


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


