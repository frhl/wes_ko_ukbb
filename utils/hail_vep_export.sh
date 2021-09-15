#!/usr/bin/env python3
import hail as hl
import argparse
import os

from gnomad.utils.vep import process_consequences
from ukb_utils import hail_init
from ko_utils import qc

def main(args):
    
    # input path
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    out_type = args.out_type
    
    # annotate with VEP
    hail_init.hail_bmrc_init('/logs/hail/hail_vep_export.log', 'GRCh38')
    dataset = qc.get_table(input_path, input_type)
    
    # clean up snpID and rsID
    dataset = qc.annotate_snpid(dataset)
    dataset = qc.annotate_rsid(dataset)
    dataset = qc.default_to_snpid_when_missing_rsid(dataset)

    # Translate to lindgren IDs
    dataset = qc.translate_sample_ids(dataset, 12788, 11867)    
    
    # recalc info
    dataset = qc.recalc_info(dataset)

    # get VEP
    result = hl.vep(dataset, "utils/configs/vep_env.json") 
    result = process_consequences(result)
    qc.export_table(result, out_prefix = out_prefix, out_type = out_type)
    
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    args = parser.parse_args()

    main(args)

