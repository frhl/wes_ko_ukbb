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

    # import phased/unphased data
    hail_init.hail_bmrc_init('logs/hail/_write_gene_intervals.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')

    # collect unique genes to investigate
    mt = io.import_table(input_path, input_type, calc_info=False)
    gene_expr = mt.consequence.vep.worst_csq_by_gene_canonical.gene_id
    genes = list(set(gene_expr.collect()))
    
    # write genes to be assessed
    idx = 1
    with open(out_prefix, "w") as outfile:
        for gene in genes:
            outfile.write(f"{idx}\t{gene}\n")
            idx += 1

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)



