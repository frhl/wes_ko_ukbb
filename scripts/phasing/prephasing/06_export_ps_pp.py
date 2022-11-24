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

    # set up
    hail_init.hail_bmrc_init_local('logs/hail/write_ps.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    mt = io.import_table(phased_path, phased_type, calc_info = False)
    mt = mt.annotate_rows(AC = mt.info.AC)
    mt = mt.annotate_rows(MAC = variants.get_mac_expr(mt))
    mt = mt.filter_entries(hl.is_defined(mt.PS_rb))
    mt = mt.select_entries(*[mt.PP, mt.GT, mt.PS_rb, mt.GT_rb])
    mt = mt.transmute_rows(rsid = variants.get_variant_expr(mt.locus, mt.alleles))
    mt = mt.select_rows(*[mt.rsid, mt.AC, mt.MAC])
    
    ht = mt.entries()
    ht.export(out_prefix + ".PP.PS.txt.gz")

    

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--phased_path', default=None, help='What is the input directory of the files?')
    parser.add_argument('--phased_type', default="vcf", help='What is the input directory of the files?')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


