#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import io
from ko_utils import ko




def main(args):
    
    input_path = args.input_path
    input_type = args.input_type
    gene = args.gene
    out_prefix = args.out_prefix
    out_type = args.out_type
    
    hail_init.hail_bmrc_init('logs/hail/gene_filter.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
    mt = io.import_table(input_path, input_type, calc_info = False)
    mt = mt.filter_rows(mt.gene_id == gene)
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='')
    parser.add_argument('--input_type', default=None, help='')
    parser.add_argument('--gene', default=None, help='')
    parser.add_argument('--out_prefix', default=None, help='')
    parser.add_argument('--out_type', default=None, help='')
    args = parser.parse_args()

    main(args)



