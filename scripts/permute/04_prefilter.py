#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import io
from ko_utils import ko

def main(args):
    
    input_path = args.input_path
    input_type = args.input_type
    csqs_category = args.csqs_category
    out_prefix = args.out_prefix
    out_type = args.out_type
    
    hail_init.hail_bmrc_init('logs/hail/prefilter.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')

    # get matrix table with variants
    mt = io.import_table(input_path, input_type)
    
    # Build variant annotation
    mt = mt.explode_rows(mt.consequence.vep.worst_csq_by_gene_canonical)
    mt = mt.annotate_rows(
        consequence_category=ko.csqs_case_builder(
                worst_csq_expr=mt.consequence.vep.worst_csq_by_gene_canonical,
                use_loftee=True)) 
    
    # 
    category = "_".join(csqs_category)
    items = csqs_category
    mt = mt.filter_rows(hl.literal(set(items)).contains(mt.consequence_category)) 
    
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='')
    parser.add_argument('--input_type', default=None, help='')
    parser.add_argument('--csqs_category', default=None, help='')
    parser.add_argument('--out_prefix', default=None, help='')
    parser.add_argument('--out_type', default=None, help='')
    args = parser.parse_args()

    main(args)



