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
    remove_withdrawn = args.remove_withdrawn
    min_info = args.min_info
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

    # get 500k calls
    calls_samples = mt.cols().s.collect()
    
    # get whole exome samples
    wes = io.import_table(in_path, in_type, calc_info=False)
    if remove_withdrawn:
        wes = samples.remove_withdrawn(wes)
    wes_samples = wes.cols().s.collect()
    
    # assume samples should be extracted
    if extract_samples:
        input_ht = hl.import_table(extract_samples, no_header=True, key='f0', delimiter=',')
        input_samples = input_ht.f0.collect()

    # get overlap
    overlap = set(calls_samples) & set(wes_samples)
    n = len(overlap)
    print(f"Note. Found {n} samples in calls and whole exomes.")
    wes = wes.filter_cols(hl.literal(overlap).contains(wes.s))
    
    # export parents that are in WES and CALLs. These are the ones
    # that we should like to evaluate downstream and thus need to exlcude.
    parents = samples.get_parents_by_fam(wes, ["TRIO"])
    mt_parents = wes.filter_cols(hl.literal(parents).contains(wes.s))
    mt_parents.cols().s.export(out_prefix + "_parents.txt")

    ht = wes.cols()
    ht.s.export(out_prefix + ".txt")
    
    # print out europeans
    wes = samples.filter_ukb_to_ancestry(wes, "eur")
    parents = samples.get_parents_by_fam(wes, ["TRIO"])
    mt_parents = wes.filter_cols(hl.literal(parents).contains(wes.s))
    mt_parents.cols().s.export(out_prefix + "_eur_parents.txt")

    ht = wes.cols()
    ht.s.export(out_prefix + "_eur.txt")
    


    

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='chromosome')
    parser.add_argument('--in_path', default=None, help='input path')
    parser.add_argument('--in_type', default=None, help='input type (vcf/mt/plink)')
    parser.add_argument('--min_info', default=None, help='imputation info when dataset is "imp"')
    parser.add_argument('--dataset', default=None, help='Either "imp" or "calls".')
    parser.add_argument('--extract_samples', default=None, help='HailTable with samples to be extracted.')
    parser.add_argument('--remove_withdrawn', default=None, action = 'store_true', help='Should UKBB withdrawn samples be discarded')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)


