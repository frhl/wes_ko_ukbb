#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import variants
from ko_utils import io
from ko_utils import ko

def main(args):
    
    in_file = args.in_file
    in_type = args.in_type
    out_prefix = args.out_prefix

    hail_init.hail_bmrc_init_local('logs/hail/get_csqs.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')

    mt = io.import_table(in_file, in_type)
    mt = io.recalc_info(mt)
    mt = mt.annotate_rows(
        MAF=variants.get_maf_expr(mt),
        MAC=variants.get_mac_expr(mt),
        consequence_category=ko.csqs_case_builder(
                worst_csq_expr=mt.consequence.vep.worst_csq_for_variant_canonical,
                use_loftee=True)) 
    ht = mt.rows()
    varid = hl.delimit([
        hl.str(ht.locus.contig), 
        hl.str(ht.locus.position), 
        ht.alleles[0], 
        ht.alleles[1]],':')
    csqs = ht.consequence.vep.worst_csq_for_variant_canonical
    ht = ht.annotate(
        varid = varid,
        csqs = csqs)
    ht = ht.select('rsid','info','varid','csqs')
    ht.flatten().export(out_prefix + ".tsv.gz")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--in_file', default=None, required = True, help='Path to input')
    parser.add_argument('--in_type', default=None, required = True, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, required = True, help='Path prefix for output dataset')
 
    args = parser.parse_args()

    main(args)

