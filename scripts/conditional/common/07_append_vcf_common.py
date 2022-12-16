#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import tables
from ko_utils import io

def main(args):

    chrom = args.chrom
    common_path = args.common_path
    common_type = args.common_type
    ko_path = args.ko_path
    ko_type = args.ko_type
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init('logs/hail/combine_ko_common.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')

    # ko only contains pseudo variants
    ko = io.import_table(ko_path, ko_type, calc_info=False)
    common = io.import_table(common_path, common_type, calc_info=False)

    # re-annotate dosage as alt alleles
    common = common.transmute_entries(DS = hl.float64(common.GT.n_alt_alleles()))

    # filter to chromosomes (Assuing common
    # variants have already been combined).
    contig = "chr" + chrom
    common = common.filter_rows(common.locus.contig == contig)

    # merge tables and export
    ko = tables.order_cols(ko, common)
    final = io.rbind_matrix_tables(ko, common)
    io.export_table(final, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='')
    parser.add_argument('--common_path', default=None, help='')
    parser.add_argument('--common_type', default=None, help='')
    parser.add_argument('--ko_path', default=None, help='')
    parser.add_argument('--ko_type', default=None, help='')
    parser.add_argument('--out_prefix', default=None, help='')
    parser.add_argument('--out_type', default=None, help='')
    args = parser.parse_args()

    main(args)
