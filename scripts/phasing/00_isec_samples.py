#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import samples

def main(args):

    # parser
    chrom = args.chrom
    dataset = args.dataset
    in_path = args.in_path
    in_type = args.in_type
    extract_samples = args.extract_samples
    min_info = args.min_info
    ancestry = args.ancestry
    out_prefix = args.out_prefix

    hail_init.hail_bmrc_init_local('logs/hail/overlapping_samples.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip

    # get imputed calls or genotyped calls
    if dataset in "imp":
        mt = genotypes.get_ukb_imputed_v3_bgen(chroms=[chrom])
        if min_info:
            ht = genotypes.get_ukb_parsed_imputed_v3_mfi(chroms=[chrom])
            mt = mt.annotate_rows(info_score=ht[(mt.locus, mt.alleles)].info)
            mt = mt.filter_rows(mt.info_score >= float(min_info))
    elif dataset in "calls":
        mt = genotypes.get_ukb_genotypes_bed(chroms=[chrom])
    else:
        raise TypeError(f"{dataset} is not 'imp' or 'calls'")

    # get calls
    calls = samples.convert_sample_ids(mt, 12788, 11867)
    calls_samples = calls.col_key.collect()
   
    print(calls.describe())

    # get whole exome samples
    wes = io.import_table(in_path, in_type, calc_info=False)
    if ancestry:
        wes = samples.filter_ukb_to_ancestry(wes, ancestry)
    wes_samples = wes.col_key.collect()

    print(wes.describe())
    
    # assume samples should be extracted
    input_ht = hl.import_table(extract_samples, no_header=True, key='f0', delimiter=',')
    input_samples = input_ht.f0.collect()
    
    
    print(calls_samples[0:5])
    print(wes_samples[0:5])
    print(input_samples[0:5])


    # get overlap
    overlap = set(calls_samples) & set(wes_samples) & set(input_samples)
    wes = wes.filter_cols(hl.literal(overlap).contains(wes.s))
    ht = wes.s
    ht.write(out_prefix + ".ht")
    ht.s.export(out_prefix + ".txt")
    

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='chromosome')
    parser.add_argument('--in_path', default=None, help='chromosome')
    parser.add_argument('--in_type', default=None, help='chromosome')
    parser.add_argument('--ancestry', default=None, help='filter to specific ancestry')
    parser.add_argument('--min_info', default=None, help='filter to specific ancestry')
    parser.add_argument('--convert_sample_id', default=None, action='store_true', help='convert to lindgren sample id')
    parser.add_argument('--dataset', default=None, help='Either "imp" or "calls".')
    parser.add_argument('--extract_samples', default=None, help='HailTable with samples to be extracted.')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


