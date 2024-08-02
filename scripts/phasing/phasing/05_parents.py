#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import variants

def main(args):
    
    children_path = args.children_path
    children_type = args.children_type
    out_prefix = args.out_prefix
    trio_path = args.trio_path

    # set up
    hail_init.hail_bmrc_init_local('logs/hail/switch_pp_subset.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulipi
    mt = io.import_table(children_path, children_type, calc_info = True) #, find_replace=(':-?nan', ':NaN'))
    
    # import output of bcftools switch errors
    ht = hl.import_table(trio_path, no_header=False)
    ht = ht.annotate(varid = hl.delimit([ht.CHR, ht.POS, ht.REF, ht.ALT],":"))
    ht = ht.key_by(**hl.parse_variant(ht.varid))
    
    # annotate with variant QC
    ht = ht.annotate(info = mt.rows().index(ht.key).info)
    ht.flatten().export(out_prefix + ".txt.gz")

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--children_path', default=None, help='What is the input directory of the files?')
    parser.add_argument('--children_type', default="vcf", help='What is the input directory of the files?')
    parser.add_argument('--trio_path', default=None, help='Path for output dataset')
    parser.add_argument('--out_prefix', default=None, help='Path for output dataset')
    args = parser.parse_args()

    main(args)


