#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import samples
from ko_utils import io

def main(args):
    
    in_file = args.in_file
    in_type = args.in_type
    out_prefix = args.out_prefix
    final_sample_list = args.final_sample_list
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/export_hets.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')

    mt = io.import_table(in_file, in_type, calc_info = True)
    
    # remove withdrawn and filter to selected samples
    if final_sample_list:
        ht_final_samples = hl.import_table(
                final_sample_list,
                no_header=True, key='f0',
                delimiter=',')
        mt = mt.filter_cols(
                hl.is_defined(ht_final_samples[mt.col_key]))

    # annotate rows
    mt = mt.annotate_rows(
        varid = hl.delimit([
            hl.str(mt.locus.contig), 
            hl.str(mt.locus.position), 
            mt.alleles[0], 
            mt.alleles[1]],':'),
        contig = hl.str(mt.locus.contig),
        position = hl.str(mt.locus.position),
        allele1 = mt.alleles[0],
        allele2 = mt.alleles[1]
    )

    # drop cols
    #mt = mt.drop("consequence")
    mt = mt.drop("filters")
    mt = mt.drop("qual")
    mt = mt.drop("rsid")
    
    # filter to hets
    mt = mt.filter_entries(mt.GT.is_het())

    # write file
    ht = mt.entries().flatten()
    ht = ht.drop("locus")
    ht = ht.drop("alleles")
    ht.export(out_prefix + ".tsv.gz")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--final_sample_list', default=None, help='Path to hail table with final samples to be included')
    parser.add_argument('--final_variant_list', default=None, help='Path to hail table with final variants to be included') 
    parser.add_argument('--in_file', default=None, required = True, help='Path to input')
    parser.add_argument('--in_type', default=None, required = True, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, required = True, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, required = True, help='either "ht" or "tsv".')
 
    args = parser.parse_args()

    main(args)

