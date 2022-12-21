#!/usr/bin/env python3

import hail as hl
import argparse
import os.path
import sys

from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import variants
from ko_utils import io

def main(args):

    chrom = args.chrom
    out_prefix = args.out_prefix
    out_type = args.out_type
    imp_path = args.imp_path
    calls_path = args.calls_path
    samples_phased_path = args.samples_phased_path
    samples_qc_path = args.samples_qc_path
    samples_extract = args.samples_extract

    reference_genome = 'GRCh38'
    hail_init.hail_bmrc_init(log='logs/hail/export_imputed.log', default_reference=reference_genome) 

    imp = hl.read_matrix_table(imp_path)
    calls = hl.read_matrix_table(calls_path)
    #ht_phased = hl.import_table(samples_phased_path, no_header = True, key = 'f0')
    #ht_qc = hl.import_table(samples_qc_path, no_header = True, key = 'f0')
    ht = hl.import_table(samples_extract, no_header = True, key = 'f0')
    calls = calls.filter_cols(hl.is_defined(ht[calls.col_key]))
    #calls = calls.filter_cols(hl.is_defined(ht_phased[calls.col_key]))
    
    cols = calls.count_cols()
    print("the cols: " + str(cols))

    # set to missing
    calls = calls.filter_entries(~hl.is_defined(calls.GT))
    calls = calls.annotate_rows(
                varid = hl.delimit(
                    [hl.str(calls.locus.contig),
                     hl.str(calls.locus.position),
                     calls.alleles[0],
                     calls.alleles[1]],
                    ':')
                ) 

    # select rows
    calls = calls.select_rows(*[calls.rsid, calls.varid])
    imp = imp.select_rows(*[imp.rsid, imp.varid])

    # export imputed
    mt = imp.union_cols(calls)
    io.export_table(mt, out_prefix, out_type)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, required=True, help='chromosome to load')
    parser.add_argument('--imp_path', default=None, required=True, help='chromosome to load')
    parser.add_argument('--calls_path', default=None, required=True, help='chromosome to load')
    parser.add_argument('--samples_phased_path', default=None, required=False, help='chromosome to load')
    parser.add_argument('--samples_qc_path', default=None, required=False, help='chromosome to load')
    parser.add_argument('--samples_extract', default=None, required=False, help='chromosome to load')
    parser.add_argument('--out_prefix', default=None, required=True, help='Path prefix for output dataset (plink format)')
    parser.add_argument('--out_type', default=None, required=True, help='Path prefix for output dataset (plink format)')
    args = parser.parse_args()

    main(args)
