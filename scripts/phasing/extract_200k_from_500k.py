#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import variants
from ukb_utils import tables

def main(args):

    # parser
    calls_path = args.calls_path
    calls_type = args.calls_type
    wes_path = args.wes_path
    wes_type = args.wes_type
    extract_samples = args.extract_samples
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/07_extract_200k.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    mt = io.import_table(calls_path, calls_type, calc_info = False) # assuming build GRCh38
    if extract_samples:
        ht_samples = hl.import_table(extract_samples, no_header=False, key='s', delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_samples[mt.col_key]))
    
    # Load WES data and exclude any overlapping variants and ensurue
    # that samples are orederd exactly the same way.
    wes = io.import_table(wes_path, wes_type, calc_info = False)
    mt = tables.order_cols(mt, wes)    
   
    #mt = mt.filter_rows(~hl.is_defined(wes.index_rows(mt.locus, mt.alleles)))

    # we can't use 'recalc_info" here, since SHAPEIT requires the AC/AN to be inputted
    # as int32 wheras the direct output of hl.agg_call_stats for AC is an array for the
    # the reference and alternate allele count. We extract the alternate allele count
    # since MAF is calculated internally in 
    mt = io.recalc_info(mt)
    
    # remove invariant sites
    invariant = (mt.info.AC == 0) | (mt.info.AC == mt.info.AN)
    mt = mt.filter_rows(~invariant)
    
    # expor table
    io.export_table(mt, out_prefix, out_type) 

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--calls_path', default=None, help='What is the input path to the file?')
    parser.add_argument('--calls_type', default=None, help='What input type?')
    parser.add_argument('--wes_path', default=None, help='What is the input path to the file?')
    parser.add_argument('--wes_type', default=None, help='What input type?')
    parser.add_argument('--extract_samples', default=None, help='HailTable with samples to be extracted.')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type out vcf/plink/mt')
    args = parser.parse_args()

    main(args)

