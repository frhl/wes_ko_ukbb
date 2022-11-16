#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ko_utils import ko
from ukb_utils import hail_init
from ukb_utils import variants


def main(args):
    
    phased_path = args.phased_path
    phased_type = args.phased_type
    out_prefix = args.out_prefix
    by = args.by
    max_maf = args.max_maf

    # set up
    hail_init.hail_bmrc_init_local('logs/hail/write_ps.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    mt = io.import_table(phased_path, phased_type, calc_info = False, find_replace=(':-?nan', ':NaN'))
   
    mt = mt.explode_rows(mt.consequence.vep[by])
    mt = mt.annotate_rows(
        MAF=variants.get_maf_expr(mt),
        MAC=variants.get_mac_expr(mt),
        consequence_category=ko.csqs_case_builder(
                worst_csq_expr=mt.consequence.vep[by],
                use_loftee=True)) 
    
    if max_maf is not None:
        mt = mt.filter_rows(mt.MAF < float(max_maf))

    mt = mt.select_entries(*[mt.PP, mt.GT])
    mt = mt.transmute_rows(rsid = variants.get_variant_expr(mt.locus, mt.alleles))
    #mt = mt.annotate_rows(AC = mt.info.AC)
    #mt = mt.annotate_rows(AF = mt.info.AF)
    mt = mt.select_rows(*[mt.rsid, mt.MAC, mt.MAF, mt.consequence_category])
    ht = mt.entries()
    ht.export(out_prefix + ".vep.PP.txt.gz")
    

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--phased_path', default=None, help='What is the input directory of the files?')
    parser.add_argument('--phased_type', default="vcf", help='What is the input directory of the files?')
    parser.add_argument('--by', default=None, required = True, help='What should be used for gene annotation')
    parser.add_argument('--max_maf', default=None, required = True, help='What should be used for gene annotation')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


