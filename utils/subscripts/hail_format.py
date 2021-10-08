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
    
    # run parser
    hail_init.hail_bmrc_init('logs/hail/hail_export.log', 'GRCh38')
    mt1 = qc.get_table(input_path=input_phased_path, input_type=input_phased_type) # 12788
    mt2 = qc.get_table(input_path=input_unphased_path, input_type=input_unphased_type) # 11867 (for singletons)

    # note: need to translate ids to combine later!
    #mt1 = qc.translate_sample_ids(mt1, 12788, 11867)
    mt1 = qc.filter_to_european(mt1, only_annotate = True)
    mt2 = qc.filter_to_european(mt2, only_annotate = True)

    ### Variant filtering/annotations
    # Using mt2 as a singleton refereence, so remove those with AC > 1
    mt2 = qc.filter_max_mac(mt2, 1)

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
        
    # annotate Gene (In the future just use vep.worst_csq_by_gene_canonical downstream..) 
    mt1 = mt1.annotate_rows(vep = mt1.vep.annotate(Gene = mt1.vep.worst_csq_by_gene_canonical.gene_id))
    mt2 = mt2.annotate_rows(vep = mt2.vep.annotate(Gene = mt2.vep.worst_csq_by_gene_canonical.gene_id))
        
    # By default add snpid id annotation
    mt1 = qc.annotate_snpid(mt1)
    mt2 = qc.annotate_snpid(mt2)
    
    # annotate mt1 with dbSNP
    mt1 = qc.annotate_rsid(mt1)
    mt1 = qc.default_to_snpid_when_missing_rsid(mt1)

    # annotate mt2 with dbSNP
    mt2 = qc.annotate_rsid(mt2)
    mt2 = qc.default_to_snpid_when_missing_rsid(mt2)
    
    # export files
    qc.export_table(mt1, out_prefix = out_prefix, out_type = 'mt')
    qc.export_table(mt2, out_prefix = out_prefix + "_singletons", out_type = 'mt')

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

