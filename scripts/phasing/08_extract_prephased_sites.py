#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init

def main(args):
    
    prephased_path = args.prephased_path
    prephased_type = args.prephased_type
    out_prefix = args.out_prefix

    # set up
    hail_init.hail_bmrc_init_local('logs/hail/merge_chunks.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    mt = io.import_table(prephased_path, prephased_type, calc_info = False)
    mt = io.recalc_info(mt)
    mt = mt.filter_entries(hl.is_defined(mt.PS))
    mt = mt.select_entries(*[mt.PS, mt.GT])
    # export mafs
    ht = mt.info.flatten()
    ht.export(out_prefix + ".maf.txt.gz")
    # export infos
    ht = mt.entries()
    ht.export(out_prefix + ".PS.txt.gz")


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--prephased_path', default=None, help='What is the input directory of the files?')
    parser.add_argument('--prephased_type', default="vcf", help='What is the input directory of the files?')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


