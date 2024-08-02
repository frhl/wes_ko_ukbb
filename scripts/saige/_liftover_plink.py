#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import variants
from ukb_utils import hail_init


def main(args):

    # parser
    liftover = args.liftover
    only_valid_contigs = args.only_valid_contigs
    
    in_prefix = args.in_prefix
    out_prefix = args.out_prefix
    in_type = args.in_type
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh37')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip

    mt = io.import_table(in_prefix, in_type, calc_info = False)

    if liftover:
        mt = variants.liftover(mt, from_build='GRCh37', to_build='GRCh38', drop_annotations=True)

    if only_valid_contigs:
        chroms = [f'chr{x}' for x in range(1,23)]
        chroms.append('chrX')
        filter_expr = hl.literal(set(chroms)).contains(mt.locus.contig)
        mt = mt.filter_rows(filter_expr)
    
    if out_type and out_prefix:
        io.export_table(mt, out_prefix, out_type)


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='chromosome')
    parser.add_argument('--liftover', default=None, action='store_true', help='perform liftover')
    parser.add_argument('--only_valid_contigs', default=None, action='store_true', help='Subset variants only normal contigs chr1..22x')
    parser.add_argument('--in_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--in_type', default=None, help='Either "mt", "vcf" or "plink"')
    parser.add_argument('--out_type', default=None, help='Either "mt", "vcf" or "plink"')
    args = parser.parse_args()

    main(args)


