
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

    maf_max = args.maf_max
    maf_min = args.maf_min
    exclude = args.exclude
    pp_cutoff = args.pp_cutoff

    # import phased/unphased data
    hail_init.hail_bmrc_init('logs/hail/knockout.log', 'GRCh38')
    mt = io.import_table(input_path, input_type, calc_info = False)
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
        mt = mt.filter_rows(~hl.literal(set(ht.varid.collect())).contains(mt.rsid))
    if pp_cutoff:
        pp_cutoff = float(pp_cutoff)
        expr_pp_cutoff = (mt.PP >= pp_cutoff) & (hl.is_defined(mt.GT))
        expr_keep = ~(hl.is_defined(mt.PP)) & (hl.is_defined(mt.GT))
        mt = mt.filter_entries((expr_pp_cutoff) | (expr_keep))

    # export result
    io.export_table(mt, out_prefix, out_type)

    # we also need the VCF
    if out_type not in "vcf":
        io.export_table(mt, out_prefix, "vcf")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--maf_min', default=None, help='Select all variants with a maf greater than the indicated values')
    parser.add_argument('--maf_max', default=None, help='Select all variants with a maf less than the indicated value')
    parser.add_argument('--exclude', default=None, help='exclude variants by rsid and/or variant id')
    parser.add_argument('--pp_cutoff', default=None, help='exclude variants by rsid and/or variant id')

    args = parser.parse_args()

    main(args)

