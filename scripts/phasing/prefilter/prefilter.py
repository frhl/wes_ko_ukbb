#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import qc
from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import samples
from ukb_utils import variants

from ko_utils.variants import filter_min_mac, filter_missing

def main(args):

    # parser
    input_path = args.input_path
    input_type = args.input_type
    extract_samples = args.extract_samples
    drop_fields = args.drop_fields
    min_mac = args.min_mac
    min_maf = args.min_maf
    missing = args.missing
    ancestry = args.ancestry
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/prefilter.py', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    mt = qc.get_table(input_path, input_type, calc_info = False) # assuming build GRCh38

    if extract_samples:
        ht_samples = hl.import_table(extract_samples, no_header=False, key='s', delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_samples[mt.col_key])) 
    if ancestry:
        mt = samples.filter_ukb_to_ancestry(mt, ancestry)
    if min_maf:
        mt = mt.filter_rows(variants.get_maf_expr(mt) >= float(min_maf))
    if min_mac:
        mt = filter_min_mac(mt, int(min_mac))
    if missing:
        before = mt.count()
        mt = filter_missing(mt, float(missing))
        after = mt.count()
        print(f"Missing filter before: {before} after: {after}")
    if missing or min_mac:
        mt = io.recalc_info(mt)
    if drop_fields:
        fields = drop_fields
        fields = fields.strip().split(",")
        mt = mt.drop(*fields)
        mt = mt.drop(mt.info)
    if out_type not in "mt":
        mt = mt.checkpoint(out_prefix + ".mt", overwrite = True)
    mt = io.recalc_info(mt)
    io.export_table(mt, out_prefix, out_type)


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='What is the input path to the file?')
    parser.add_argument('--input_type', default=None, help='What input type?')
    parser.add_argument('--ancestry', default=None, help='filter to specific ancestry')
    parser.add_argument('--extract_samples', default=None, help='HailTable with samples to be extracted.')
    parser.add_argument('--min_mac', default=None, help='Filter to MAC >= value')
    parser.add_argument('--min_maf', default=None, help='Filter to MAC >= value')
    parser.add_argument('--missing', default=None, help='Filter to variants to have le value in genotype missingness')
    parser.add_argument('--drop_fields', default=None, help='Filter to variants to have le value in genotype missingness')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type out vcf/plink/mt')
    args = parser.parse_args()

    main(args)

