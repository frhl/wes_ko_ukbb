#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import variants

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
    
    mt = mt.transmute_rows(rsid = variants.get_variant_expr(mt.locus, mt.alleles))
    mt = mt.annotate_rows(AC = mt.info.AC)
    mt = mt.annotate_rows(AN = mt.info.AN)
    mt = mt.select_rows(*[mt.rsid, mt.AC, mt.AN])
    
    # get entries
    ht = mt.entries()
    ht = ht.drop(*[ht.locus, ht.alleles])
    ht.export(out_prefix + ".PS.txt.gz")

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--prephased_path', default=None, help='What is the input directory of the files?')
    parser.add_argument('--prephased_type', default="vcf", help='What is the input directory of the files?')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


