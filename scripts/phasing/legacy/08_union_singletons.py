#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import tables
from ukb_utils import variants
from ko_utils import io
from ko_utils import variants as ko_utils_variants


def main(args):
    input_phased = args.input_phased
    input_phased_type = args.input_phased_type
    input_singletons = args.input_singletons
    input_singletons_type = args.input_singletons_type
    exlucde_calls = args.exclude_calls
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/union_singletons.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')

    # load phased data
    phased = io.import_table(
            input_path=input_phased,
            input_type=input_phased_type,
            force_bgz = False
            )
    
    # all singletons that are not read-backed
    expr_singleton = variants.get_MAC_expr(phased) == 1

    # create a new matrixtable with singletons?
    # set all singletons without a PS tag to unphased?
    # no need to fold back in signletons externally.

    singletons = io.import_table(
            input_path=input_singletons,
            input_type=input_singletons_type)

    singletons = singletons.annotate_rows(wes = 1)
    singletons = ko_utils_variants.filter_max_mac(singletons, 1)
    singletons = tables.order_cols(singletons, phased)
    mt = io.rbind_matrix_tables(singletons, phased)
    if exclude_calls:
        mt = mt.filter_rows(mt.wes == 1)
    io.export_table(mt, out_prefix, out_type)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_phased', default=None, help='Path to input')
    parser.add_argument('--input_phased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_singletons', default=None, help='Path to input that contains singletons')
    parser.add_argument('--input_singletons_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_annotation_path', default=None, help='Input for annotation path')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--exclude_calls', default=None, action='store_true', help='exclude genotype array calls that were used for phasing')
    parser.add_argument('--chrom', default=None, help='Chromosome to be used')
 
    args = parser.parse_args()

    main(args)

