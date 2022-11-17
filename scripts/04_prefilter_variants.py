#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import samples
from ukb_utils import variants
from ko_utils import io
from ko_utils import ko

def main(args):

    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    out_type = args.out_type

    sex = args.sex
    maf_max = args.maf_max
    maf_min = args.maf_min
    exclude = args.exclude
    use_loftee = args.use_loftee
    pp_cutoff = args.pp_cutoff

    # import phased/unphased data
    hail_init.hail_bmrc_init('logs/hail/knockout.log', 'GRCh38')
    mt = io.import_table(input_path, input_type, calc_info = True)

    # some filtering 
    if sex not in 'both':
        mt = samples.filter_to_sex(mt, sex)
    if maf_min:
        maf_min = float(maf_min)
        if maf_min > 0:
            maf_min_expr = variants.get_maf_expr(mt)
            mt = mt.filter_rows(maf_min_expr >= maf_min)
    if maf_max:
        maf_max = float(maf_max)
        if maf_max < 1:
            maf_max_expr = variants.get_maf_expr(mt)
            mt = mt.filter_rows(maf_max_expr <= maf_max)
    if exclude:
        ht = hl.import_table(exclude, impute=True).key_by('varid')
        mt = mt.filter_rows(~hl.literal(set(ht.varid.collect())).contains(mt.varid))
    if pp_cutoff:
        pp_cutoff = float(pp_cutoff)
        mt = mt.filter_entries(mt.PP <= pp_cutoff)

    # explode by rows
    mt = mt.explode_rows(mt.consequence.vep.worst_csq_by_gene_canonical)
    mt = mt.annotate_rows(
        consequence_category=ko.csqs_case_builder(
                worst_csq_expr=mt.consequence.vep.worst_csq_by_gene_canonical,
                use_loftee=use_loftee))

    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--sex', default='both', help='Filter to sex (males or females)')
    parser.add_argument('--maf_min', default=None, help='Select all variants with a maf greater than the indicated values')
    parser.add_argument('--maf_max', default=None, help='Select all variants with a maf less than the indicated value')
    parser.add_argument('--exclude', default=None, help='exclude variants by rsid and/or variant id')
    parser.add_argument('--pp_cutoff', default=None, help='exclude variants by rsid and/or variant id')
    parser.add_argument('--use_loftee', default=False, action='store_true', help='use LOFTEE to distinghiush between high confidence PTVs')

    args = parser.parse_args()

    main(args)



