#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import io

def main(args):

    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init('logs/hail/knockout.log', 'GRCh38')
    mt = io.import_table(input_path, input_type, calc_info = True) 
    
    # remove invariant sites
    mt = mt.filter_rows(mt.info.AC>0)
    
    # annotate with variant ID
    mt = mt.transmute_rows(
                varid = hl.delimit(
                    [hl.str(mt.locus.contig),
                     hl.str(mt.locus.position),
                     mt.alleles[0],
                     mt.alleles[1]],
                    ':')
                )

    # annotate ID
    mt = mt.transmute_rows(ID=mt.varid)
   

    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')

    args = parser.parse_args()

    main(args)



