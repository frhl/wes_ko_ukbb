#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import qc
from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import samples

from ko_utils.variants import filter_min_mac, filter_missing

def main(args):

    # parser
    input_path = args.input_path
    input_type = args.input_type
    extract_samples = args.extract_samples
    exclude_trio_parents = args.exclude_trio_parents
    export_parents = args.export_parents
    drop_entry_fields = args.drop_entry_fields
    min_mac = args.min_mac
    missing = args.missing
    ancestry = args.ancestry
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/01_prefilter_wes.py', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    mt = qc.get_table(input_path, input_type, calc_info = False) # assuming build GRCh38

    if split_parents:
        mt = mt.checkpoint(out_prefix + ".mt")
        io.export_table(mt, out_prefix, out_type)
        pids = samples.get_parents_by_fam(mt, ["TRIO"])
        if export_parents:
            mt_parents = mt.filter_cols(hl.literal(pids).contains(mt.s))
            io.export_table(mt_parents, out_prefix + "_parents", out_type)
        mt = mt.filter_cols(~hl.literal(pids).contains(mt.s))

    if exclude_trio_parents:
        io.export_table(mt, out_prefix + "_no_parents", out_type) 
    else:
        io.export_table(mt, out_prefix, out_type) 


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='What is the input path to the file?')
    parser.add_argument('--input_type', default=None, help='What input type?')
    parser.add_argument('--split_parents', default=None, help='Exclude parents of trio relationships')
    parser.add_argument('--export_parents', default=None, action='store_true', help='Exclude parents of trio relationships')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type out vcf/plink/mt')
    args = parser.parse_args()

    main(args)

