#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import variants

def main(args):

    # parser
    input_path = args.input_path
    input_type = args.input_type
    extract_samples = args.extract_samples
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/07_extract_200k.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    mt = io.import_table(input_path, input_type, calc_info = False) # assuming build GRCh38
    if extract_samples:
        ht_samples = hl.import_table(extract_samples, no_header=False, key='s', delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_samples[mt.col_key]))
    
    # we can't use 'recalc_info" here, since SHAPEIT requires the AC/AN to be inputted
    # as int32 wheras the direct output of hl.agg_call_stats for AC is an array for the
    # the reference and alternate allele count. In this case, we will assume that the AC
    # is actually the minor allele count.
    mt = mt.drop(mt.info)
    mt = mt.annotate_rows(info=hl.struct())
    mt = mt.annotate_rows(info=mt.info.annotate(AC=variants.get_mac_expr(mt)))
    mt = mt.annotate_rows(info=mt.info.annotate(AN=hl.agg.call_stats(mt.GT, mt.alleles).AN))
    io.export_table(mt, out_prefix, out_type) 

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='What is the input path to the file?')
    parser.add_argument('--input_type', default=None, help='What input type?')
    parser.add_argument('--extract_samples', default=None, help='HailTable with samples to be extracted.')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type out vcf/plink/mt')
    args = parser.parse_args()

    main(args)

