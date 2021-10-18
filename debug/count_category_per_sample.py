#!/usr/bin/env python3

import hail as hl
import argparse
import pandas
import os

from gnomad.utils.vep import process_consequences
from ukb_utils import hail_init
from ukb_utils import genotypes
from ko_utils import qc
from ko_utils import analysis

def main(args):
    
    # parser
    input_phased_path = args.input_phased_path
    input_phased_type = args.input_phased_type
    input_unphased_path = args.input_unphased_path
    input_unphased_type = args.input_unphased_type
    out_prefix = args.out_prefix
    out_type   = args.out_type
    vep_path   = args.vep_path
    chrom      = args.chrom

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    
    # get tables
    mt1 = qc.get_table(input_path=input_phased_path, input_type=input_phased_type) # 12788
    mt2 = qc.get_table(input_path=input_unphased_path, input_type=input_unphased_type) # 11867 (for singletons)

    # Run VEP + gnomAD variant annotations
    mt1 = process_consequences(hl.vep(mt1, "utils/configs/vep_env.json"))
    mt2 = process_consequences(hl.vep(mt2, "utils/configs/vep_env.json"))
    
    # Annotate with REVEL+CADD scores
    mt1 = analysis.annotate_dbnsfp(mt1, vep_path)
    mt2 = analysis.annotate_dbnsfp(mt2, vep_path)
   
    # Annotate for each consequence
    mt1 = mt1.explode_rows(mt1.vep.worst_csq_by_gene_canonical)
    mt2 = mt2.explode_rows(mt2.vep.worst_csq_by_gene_canonical)

    # annotate consequnece categories 
    mt1 = analysis.variant_csqs_category_builder(mt1)
    mt2 = analysis.variant_csqs_category_builder(mt2)
    
    # filter to PTVs
    mt1 = mt1.filter_rows(mt1.vep.consequence_category == 'ptv')
    mt2 = mt2.filter_rows(mt2.vep.consequence_category == 'ptv')

    # counts
    mt1 = mt1.annotate_entries(homozygous = mt1.GT.is_hom_var())
    mt2 = mt2.annotate_entries(homozygous = mt2.GT.is_hom_var())

    mt1 = mt1.annotate_entries(heterozygous = mt1.GT.is_het_ref())
    mt2 = mt2.annotate_entries(heterozygous = mt2.GT.is_het_ref())
    
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_phased_path', default=None, help='Path to input')
    parser.add_argument('--input_phased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_unphased_path', default=None, help='Path to input that contains singletons')
    parser.add_argument('--input_unphased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--chrom', default=None, help='Chromosome to be used')
    parser.add_argument('--vep_path', default=None, help='path to a .vcf file containing annotated entries by locus and alleles')   
 
    args = parser.parse_args()

    main(args)

