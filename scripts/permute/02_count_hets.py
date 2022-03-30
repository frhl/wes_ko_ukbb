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

    # run parser
    hail_init.hail_bmrc_init('logs/hail/aggr_hets.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') 
    mt = io.import_table(input_path, input_type) 
    expr = mt.consequence.vep.worst_csq_by_gene_canonical.gene_id
    # aggregate counts
    mt = (mt.group_rows_by(expr)
                    .aggregate(
                        unphased_het=hl.agg.count_where((mt.GT.is_het_ref()) & (~mt.GT.phased)),
                        phased_het=hl.agg.count_where((mt.GT.is_het_ref()) & (mt.GT.phased)),
                        hom_alt_n=hl.agg.count_where(mt.GT.is_hom_var())
                        )
               )
    # calculate likelihood of being a true knockouts
    pTKO = (hl.case()
        .when((mt.hom_alt_n > 0), 1)
        .when((mt.phased_het > 1), 1 - 2*(1/2) ** mt.phased_het)
        .default(0))
    mt = mt.annotate_entries(pTKO = pTKO) 
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    args = parser.parse_args()

    main(args)



