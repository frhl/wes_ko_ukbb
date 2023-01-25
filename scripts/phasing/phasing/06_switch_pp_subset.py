#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import variants

def main(args):
    
    phased_path = args.phased_path
    phased_type = args.phased_type
    out_prefix = args.out_prefix
    pp_cutoff = args.pp_cutoff

    # set up
    hail_init.hail_bmrc_init_local('logs/hail/switch_pp_subset.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulipi
    # note, we don't recalculate INFO after including parents
    # as this will make singletons into doubletons/tripletons.
    mt = io.import_table(phased_path, phased_type, calc_info = False, find_replace=(':-?nan', ':NaN'))
    # note that PPs are only af MAF < 0.1%
    mt = mt.select_entries(*[mt.PP, mt.GT])
    # subset entries that are either not defined or wihtin cutoff 
    pp_cutoff = float(pp_cutoff)
    expr_pp_cutoff = (mt.PP >= pp_cutoff) & (hl.is_defined(mt.GT))
    expr_keep = ~(hl.is_defined(mt.PP)) & (hl.is_defined(mt.GT))
    mt = mt.filter_entries((expr_pp_cutoff) | (expr_keep)) 
    mt = mt.select_rows(mt.rsid)
    io.export_table(mt, out_prefix, out_type = "vcf")


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--phased_path', default=None, help='What is the input directory of the files?')
    parser.add_argument('--phased_type', default="vcf", help='What is the input directory of the files?')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--pp_cutoff', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


