#!/usr/bin/env python3

import hail as hl
import argparse
import pandas
import os

from gnomad.utils.vep import process_consequences
from ukb_utils import hail_init
from ko_utils import io


def main(args):

    # parser
    json_path = args.json_path
    input_path = args.input_path
    out_prefix = args.out_prefix

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')  # from zulip
    if "bim" in input_path:
        ht = hl.import_table(input_path, no_header=True)
        ht = ht.rename({'f0': 'contig', 'f1': 'rsid', 'f2': 'cm_position',
                        'f3': 'position', 'f4': 'allele1', 'f5': 'allele2'})
        ht = ht.key_by(locus=hl.locus(ht.contig, hl.int32(
            ht.position), reference_genome='GRCh38'), alleles=[ht.allele1, ht.allele2])
    elif "mt" in input_path:
        mt = hl.read_matrix_table(input_path)
        ht = mt.rows()
    elif "vcf" in input_path:
        mt = hl.import_vcf(input_path, reference_genome='GRCh38', force_bgz=True)
        ht = mt.rows().select()
    else:
        ht = hl.read_table(input_path)

    ht = hl.vep(ht, json_path, tolerate_parse_error=False)
    ht.export(out_prefix + ".txt.gz")
    ht.write(out_prefix + ".ht", overwrite=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--json_path', default=None, help='Path to json')
    parser.add_argument('--tolerate_parse_error', default=False, action='store_true')
    parser.add_argument('--out_prefix', default=None,
                        help='Path prefix for output dataset')

    args = parser.parse_args()

    main(args)

