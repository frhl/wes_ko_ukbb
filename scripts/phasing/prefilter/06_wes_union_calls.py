#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ko_utils import ko
from ukb_utils import hail_init
from ukb_utils import tables

def main(args):

    # parser
    chrom = args.chrom
    input_calls_path = args.input_calls_path
    input_calls_type = args.input_calls_type
    input_wes_path = args.input_wes_path
    input_wes_type = args.input_wes_type
    out_calls_prefix = args.out_calls_prefix
    out_calls_type = args.out_calls_type
    extract_samples_calls = args.extract_samples_calls
    unphase = args.unphase

    hail_init.hail_bmrc_init_local('logs/hail/wes_union_calls_bcf.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip

    # load calls and whole exomes    
    mt1 = io.import_table(input_calls_path, input_calls_type, calc_info = False)
    mt2 = io.import_table(input_wes_path, input_wes_type, calc_info = False) # assuming build GRCh38
    
    # ensure that overlapping columns are retained
    if extract_samples_calls:
        ht_samples_calls = hl.import_table(extract_samples_calls, no_header=False, key='s', delimiter=',')
        mt1 = mt1.filter_cols(hl.is_defined(ht_samples_calls[mt1.col_key]))

    #mt1 = mt1.filter_cols(~hl.is_defined(mt2.index_cols(mt1.s)))
    
    # in order to combine with bcftools, we need to order by samples
    mt1 = tables.order_cols(mt1, mt2)

    # only keep keys and GT for each MatrixTable
    mt1 = mt1.drop(*set(list(mt1.col)) - set(list(mt1.col_key)))
    mt1 = mt1.drop(*set(list(mt1.row)) - set(list(mt1.row_key)))
    mt1 = mt1.select_entries(mt1.GT)

    mt2 = mt2.drop(*set(list(mt2.col)) - set(list(mt2.col_key)))
    mt2 = mt2.drop(*set(list(mt2.row)) - set(list(mt2.row_key)))
    mt2 = mt2.select_entries(mt2.GT)  

    # Remove any variants from mt that are already in mt2
    mt1 = mt1.filter_rows(~hl.is_defined(mt2.index_rows(mt1.locus, mt1.alleles)))
    
    # unphase entries
    if unphase:
        mt1 = mt1.transmute_entries(GT = ko.unphase(mt1.GT))

    # recalculate INFO
    mt1 = io.recalc_info(mt1)
    io.export_table(mt1, out_calls_prefix, out_calls_type)

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='chromosome')
    parser.add_argument('--input_wes_path', default=None, help='What is the input path to the file?')
    parser.add_argument('--input_calls_path', default=None, help='What is the input path to the file?')
    parser.add_argument('--input_wes_type', default=None, help='What input type?')
    parser.add_argument('--input_calls_type', default=None, help='What input type?')
    parser.add_argument('--out_calls_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_calls_type', default=None, help='Type out vcf/plink/mt')
    parser.add_argument('--unphase', default=None, action='store_true', help='Unphase genotypes')
    parser.add_argument('--extract_samples_calls', default=None, help='extract samples from CALLs')
    args = parser.parse_args()

    main(args)


