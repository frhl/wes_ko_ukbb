#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import tables
from ko_utils import io
from ko_utils import variants as ko_utils_variants


def main(args):
    input_phased = args.input_phased
    input_phased_type = args.input_phased_type
    input_unphased = args.input_unphased
    input_unphased_type = args.input_unphased_type
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/union_mts.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')

    phased = io.import_table(
            input_path=input_phased,
            input_type=input_phased_type)

    unphased = io.import_table(
            input_path=input_unphased,
            input_type=input_unphased_type)

    unphased = unphased.annotate_rows(wes = 1)
    unphased = ko_utils_variants.filter_max_mac(unphased, 1)
    unphased = tables.order_cols(unphased, phased)
    mt = io.rbind_matrix_tables(unphased, phased)
    io.export_table(mt, out_prefix, out_type)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_phased', default=None, help='Path to input')
    parser.add_argument('--input_phased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_unphased', default=None, help='Path to input that contains singletons')
    parser.add_argument('--input_unphased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_annotation_path', default=None, help='Input for annotation path')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--chrom', default=None, help='Chromosome to be used')
 
    args = parser.parse_args()

    main(args)

