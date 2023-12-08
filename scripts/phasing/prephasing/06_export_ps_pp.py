#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import variants

def main(args):

    phased_path = args.phased_path
    phased_type = args.phased_type
    out_prefix = args.out_prefix
    sample_size = args.sample_size
    seed = args.seed
    sample_file = args.sample_file

    # set up
    hail_init.hail_bmrc_init_local('logs/hail/write_ps.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')  # from zulip
    mt = io.import_table(phased_path, phased_type, calc_info=False)

    # Subset to specific samples if a sample file is provided
    if sample_file:
        with open(sample_file, 'r') as f:
            sample_ids = [line.strip() for line in f]
        mt = mt.filter_cols(hl.literal(sample_ids).contains(mt.s))

    # Then, subset to a random sample if sample size is specified
    if sample_size:
        mt = mt.sample_rows(p=sample_size/mt.count_rows(), seed=seed)

    # grab the allele counts from the WES200k
    mt = mt.annotate_rows(AC=mt.info.AC[0])
    mt = mt.annotate_rows(AN=mt.info.AN)
    mt = mt.annotate_rows(AF=mt.info.AF)
    mt = mt.annotate_rows(MAC=hl.min(mt.AC, mt.AN-mt.AC))
   
    # we now filter to defined phased sets in the data
    mt = mt.filter_entries(hl.is_defined(mt.PS_rb))
    mt = mt.select_entries(*[mt.PP, mt.GT, mt.PS_rb, mt.GT_rb])
    
    # re-annotate variants 
    mt = mt.transmute_rows(rsid=variants.get_variant_expr(mt.locus, mt.alleles))
   
    # keep wes200k AC/AN and also the AC for the read-backed phased data
    mt = mt.select_rows(*[mt.rsid, mt.AC, mt.AN, mt.MAC])

    ht = mt.entries()
    ht.export(out_prefix + ".PP.PS.txt.gz")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--phased_path', default=None, help='What is the input directory of the files?')
    parser.add_argument('--phased_type', default="vcf", help='What is the input directory of the files?')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--sample_size', type=int, default=None, help='Number of samples to subset (e.g., 10000)')
    parser.add_argument('--seed', type=int, default=42, help='Seed for reproducibility')
    parser.add_argument('--sample_file', default=None, help='File containing sample IDs for subsetting')
    args = parser.parse_args()

    main(args)

