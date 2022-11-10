#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init

def main(args):
    
    phased_path = args.phased_path
    phased_type = args.phased_type
    out_prefix = args.out_prefix
    max_mac = args.max_mac

    # set up
    hail_init.hail_bmrc_init_local('logs/hail/merge_chunks.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    phased = io.import_table(phased_path, phased_type, calc_info = False)
    phased = io.recalc_info(phased)
    phased = phased.annotate_rows(min_mac = hl.min(phased.info.AC, phased.info.AN - phased.info.AC))
    phased = phased.filter_rows(phased.min_mac <= int(max_mac))
    phased = phased.select_entries(phased.PP)
    phased = phased.select_rows(phased.min_mac)
    ht = phased.make_table()
    ht = ht.flatten()
    ht.export(out_prefix + ".txt.gz")


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--phased_path', default=None, help='What is the input directory of the files?')
    parser.add_argument('--phased_type', default="vcf", help='What is the input directory of the files?')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--max_mac', default=50, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


