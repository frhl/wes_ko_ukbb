#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import samples

from ko_utils.variants import filter_min_mac, filter_missing

def main(args):

    # parser
    input_path = args.input_path
    input_type = args.input_type
    parents_path = args.parents_path
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/01_prefilter_wes.py', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    mt = io.import_table(input_path, input_type, calc_info = False) # assuming build GRCh38
    #pids = samples.get_parents_by_fam(mt, ["TRIO"])
    pids = hl.import_table(parents_path, no_header=False, key='s', 
            delimiter=',', types={'s': hl.tstr}).s.collect()   
    mt_parents = mt.filter_cols(hl.literal(pids).contains(mt.s))
    io.export_table(mt_parents, out_prefix + "_parents", out_type)
    mt = mt.filter_cols(~hl.literal(pids).contains(mt.s))
    io.export_table(mt, out_prefix + "_no_parents", out_type) 


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='What is the input path to the file?')
    parser.add_argument('--input_type', default=None, help='What input type?')
    parser.add_argument('--parents_path', default=None, help='What input type?')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type out vcf/plink/mt')
    args = parser.parse_args()

    main(args)

