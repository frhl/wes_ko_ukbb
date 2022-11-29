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
    hl._set_flags(no_whole_stage_codegen='1') # from zulipi
    # note, we don't recalculate INFO after including parents
    mt = io.import_table(phased_path, phased_type, calc_info = True, find_replace=(':-?nan', ':NaN'))
    # note that PPs are only af MAF < 0.1%
    mt = mt.select_entries(*[mt.PP, mt.GT])
    # we are only intersted in heterozygous sites with PP field defined
    mt = mt.filter_entries(hl.is_defined(mt.PP))
    mt = mt.filter_entries(mt.GT.is_het())
    mt = mt.transmute_rows(rsid = variants.get_variant_expr(mt.locus, mt.alleles))
    mt = mt.annotate_rows(MAC = variants.get_mac_expr(mt))
    mt = mt.annotate_rows(AC = mt.info.AC)
    mt = mt.annotate_rows(AN = mt.info.AN)
    mt = mt.select_rows(*[mt.rsid, mt.MAC, mt.AC, mt.AN])
    ht = mt.entries()
    ht.export(out_prefix + ".qc.conf.txt.gz")
    

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--phased_path', default=None, help='What is the input directory of the files?')
    parser.add_argument('--phased_type', default="vcf", help='What is the input directory of the files?')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


