#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import tables
from ko_utils import io

def main(args):
    
    ko_path = args.ko_path
    ko_type = args.ko_type
    markers_path = args.markers_path
    markers_type = args.markers_type
    markers = args.markers
    out_prefix = args.out_prefix
    out_type = args.out_type
    
    hail_init.hail_bmrc_init('logs/hail/spa_conditional.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
    markers = markers.split(",") # format chr12:4214:A:T
    mk = io.import_table(markers_path, markers_type, calc_info = False)
    mk = mk.annotate_rows(
            marker = hl.delimit([
                hl.str(mk.locus.contig), 
                hl.str(mk.locus.position),
                hl.str(mk.alleles[0]),
                hl.str(mk.alleles[1])
                ],':'))
    mk = mk.filter_rows(hl.literal(set(markers)).contains(mk.marker))
    # combine the two tables and export
    ko = io.import_table(ko_path, ko_type, calc_info = False)
    ko = tables.order_cols(ko, mk)
    mt = io.rbind_matrix_tables(ko, mk)
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='')
    parser.add_argument('--input_type', default=None, help='')
    parser.add_argument('--markers_path', default=None, help='')
    parser.add_argument('--markers_type', default=None, help='')
    parser.add_argument('--markers', default=None, help='')
    parser.add_argument('--out_prefix', default=None, help='')
    parser.add_argument('--out_type', default=None, help='')
    args = parser.parse_args()

    main(args)



