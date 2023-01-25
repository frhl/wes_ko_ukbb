#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import io


class SplitArgs(argparse.Action):
    r"""Method for splitting input csv into a list"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))

def main(args):
    
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    out_type = args.out_type 
    csqs_category = args.csqs_category 

    # run parser
    hail_init.hail_bmrc_init('logs/hail/aggr_hets.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') 
    mt = io.import_table(input_path, input_type) 
    
    # Build variant annotation
    mt = mt.filter_rows(hl.literal(set(csqs_category)).contains(mt.consequence_category))
    
    # aggregate counts`
    expr = mt.consequence.vep.worst_csq_by_gene_canonical.gene_id
    mt = (mt.group_rows_by(expr)
                    .aggregate(
                        het_n=hl.agg.count_where((mt.GT.is_het_ref()) & (mt.GT.phased)),
                        hom_alt_n=hl.agg.count_where(mt.GT.is_hom_var())
                        )
               )
    
    # What is the probability that they are all on the same phase? We don't need to distinghuish
    # between mac1 an mac+1 PTVs in this case. 
    # 1 - p(all on phase 1) - p(all on phase 2) = 1 - 2*p(all on phase 1) = 1 - 2*(1/2)^k
    
    # calculate likelihood of being a true knockouts, these probabilities
    # will be used later to draw from when permuting the phase.
    pTKO = (hl.case()
        .when((mt.hom_alt_n > 0), 1)
        .when((mt.het_n > 1), 1 - 2*(1/2) ** mt.het_n)
        .default(0))
    mt = mt.annotate_entries(pTKO = pTKO) 
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='What categories should be subsetted to?')
    
    args = parser.parse_args()


    main(args)



