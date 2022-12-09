#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ko_utils import ko
from ukb_utils import hail_init
from ukb_utils import samples
from ukb_utils import variants

from ko_utils.variants import filter_min_mac, filter_missing

def main(args):

    # parser
    input_path = args.input_path
    input_type = args.input_type
    parents_path = args.parents_path
    no_metadata = args.no_metadata
    out_prefix = args.out_prefix
    out_type = args.out_type
    unphase = args.unphase

    hail_init.hail_bmrc_init_local('logs/hail/01_prefilter_wes.py', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    mt = io.import_table(input_path, input_type, calc_info = False) # assuming build GRCh38
    if unphase:
        mt = mt.transmute_entries(GT = ko.unphase(mt.GT))
    pids = hl.import_table(parents_path, no_header=False, key='s', 
            delimiter=',', types={'s': hl.tstr}).s.collect()   
    
    # filter to subset of samples that are parents
    mt_parents = mt.filter_cols(hl.literal(pids).contains(mt.s))
    #io.export_table(mt_parents, out_prefix + "_parents", "mt", auto_metadata=no_metadata)
    io.export_table(mt_parents, out_prefix + "_parents", out_type, auto_metadata=no_metadata)
    
    # filter to samples without parents
    #mt = mt.filter_cols(~hl.literal(pids).contains(mt.s))
    #mt = mt.filter_rows(~variants.get_invariant_expr(mt))
    #mt = io.recalc_info(mt)
    #io.export_table(mt, out_prefix + "_no_parents", "mt", auto_metadata=no_metadata) 
    #io.export_table(mt, out_prefix + "_no_parents", out_type, auto_metadata=no_metadata) 


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='What is the input path to the file?')
    parser.add_argument('--input_type', default=None, help='What input type?')
    parser.add_argument('--parents_path', default=None, help='What input type?')
    parser.add_argument('--no_metadata', default=True, action='store_false', help='What input type?')
    parser.add_argument('--unphase', default=True, action='store_true', help='What input type?')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type out vcf/plink/mt')
    args = parser.parse_args()

    main(args)

