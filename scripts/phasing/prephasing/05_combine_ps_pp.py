#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init

def main(args):
    
    phased_path = args.phased_path
    phased_type = args.phased_type
    ref_path = args.ref_path
    ref_type = args.ref_type    
    out_prefix = args.out_prefix
    out_type = args.out_type

    # set up
    hail_init.hail_bmrc_init_local('logs/hail/merge_chunks.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    # original 200k phased data
    mt1 = hl.import_vcf(
        phased_path,
        force_bgz=True,
        skip_invalid_loci=True,
        array_elements_required=False,
        find_replace=(':-nan', ':NaN'))
    
    # readback phased data
    mt2 = hl.import_vcf(
        ref_path,
        force_bgz=True,
        skip_invalid_loci=True,
        array_elements_required=False,
        find_replace=('nan', 'NaN'))
    
    # we only care about read-backed entries in 
    # which we have a phased set (i.e. phased genotypes)
    mt2 = mt2.filter_entries(hl.is_defined(mt2.PS))
    
    # transfer the read-backed phased entries to
    # the original phased data
    mt1 = mt1.annotate_entries(
        GT_rb = mt2[mt1.row_key, mt1.col_key].GT,
        PS_rb = mt2[mt1.row_key, mt1.col_key].PS
    )   

    # we also keep the allele count from the read-backed
    # phased data. Just in case, we want to filter downstream
    mt1.transmute_rows(
        info = mt1.info.annotate(
            AC_rb = mt2.rows()[mt1.row_key].info.AC,
            AN_rb = mt2.rows()[mt1.row_key].info.AN
        )
    )

    io.export_table(mt1, out_prefix, out_type)


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--phased_path', default=None, help='What is the input directory of the files?')
    parser.add_argument('--phased_type', default="vcf", help='What is the input directory of the files?')
    parser.add_argument('--ref_path', default=None, help='What is the input directory of the files?')
    parser.add_argument('--ref_type', default="vcf", help='What is the input directory of the files?')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


