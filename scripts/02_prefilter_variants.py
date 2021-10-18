#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import qc

def main(args):
    
    # parser
    input_phased_path = args.input_phased_path
    input_phased_type = args.input_phased_type
    input_unphased_path = args.input_unphased_path
    input_unphased_type = args.input_unphased_type
    out_prefix = args.out_prefix
    out_type   = args.out_type
    vep_path   = args.vep_path
    
    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    
    # get tables
    mt1 = qc.get_table(input_path=input_phased_path, input_type=input_phased_type)
    mt2 = qc.get_table(input_path=input_unphased_path, input_type=input_unphased_type) 

    # perform variant QC    
    mt1 = hl.variant_qc(mt1, name='variant_qc')
    mt2 = hl.variant_qc(mt2, name='variant_qc')
    mt1 = mt1.filter_rows((mt1.variant_qc.AF[0] > 0.0) & (mt1.variant_qc.AF[0] < 1.0))
    mt2 = mt2.filter_rows((mt2.variant_qc.AF[0] > 0.0) & (mt2.variant_qc.AF[0] < 1.0))

    # write out to file
    ht_rows_filter1 = mt1.rows()
    ht_rows_filter2 = mt2.rows()
    ht_rows_filter1.select().write(out_filter_prefix + "_phased"), overwrite=True)
    ht_rows_filter2.select().write(out_filter_prefix + "_unphased"), overwrite=True)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_phased_path', default=None, help='Path to input')
    parser.add_argument('--input_phased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_unphased_path', default=None, help='Path to input that contains singletons')
    parser.add_argument('--input_unphased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
 
    args = parser.parse_args()

    main(args)

