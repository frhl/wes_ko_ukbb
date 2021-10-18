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
    initial_variant_list = args.initial_variant_list

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    mt = qc.get_table(input_path=input_path, input_type=input_type)

    # load variants
    variants_to_filter = hl.read_table(initial_variant_list)
    
    # filter variants and get sample QC
    mt = mt.filter_rows(hl.is_defined(variants_to_filter[mt.row_key]))
    mt = mt.annotate_cols(gq = hl.agg.stats(mt.GQ), dp = hl.agg.stats(mt.DP))
    mt = hl.sample_qc(mt, name='sample_qc')

    mt.cols().select('sample_qc', 'gq', 'dp').flatten().export(output=out_prefix + '_initial_sample_qc.tsv.bgz')

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--initial_variant_list', default=None, help='List of initial QCed variants')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
 
    args = parser.parse_args()

    main(args)

