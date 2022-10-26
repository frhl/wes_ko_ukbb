#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import qc
from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import tables
from ukb_utils import genotypes
from ukb_utils import variants
from ukb_utils import samples
from ko_utils.samples import filter_to_females
from ko_utils.variants import filter_min_mac, filter_missing

def main(args):

    # parser
    chrom = args.chrom
    input_path = args.input_path
    input_type = args.input_type
    chrom = args.chrom
    extract_samples = args.extract_samples
    exclude_trio_parents = args.exclude_trio_parents
    export_parents = args.export_parents
    min_mac = args.min_mac
    missing = args.missing
    ancestry = args.ancestry
    checkpoint = args.checkpoint
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/01_geno_gen.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    print("Current Chrom:" + str(chrom))

    if input_path and input_type:

        mt2 = qc.get_table(input_path, input_type) # assuming build GRCh38
        
        # only keep keys and GT for each MatrixTable
        mt = mt.drop(*set(list(mt.col)) - set(list(mt.col_key)))
        mt = mt.drop(*set(list(mt.row)) - set(list(mt.row_key)))
        mt = mt.select_entries(mt.GT)

        mt2 = mt2.drop(*set(list(mt2.col)) - set(list(mt2.col_key)))
        mt2 = mt2.drop(*set(list(mt2.row)) - set(list(mt2.row_key)))
        mt2 = mt2.select_entries(mt2.GT)  

        # subset to intersecting samples
        mt_sids = mt.s.collect()
        mt2_sids = mt2.s.collect()
        overlap = list(set(mt_sids) & set(mt2_sids))
        
        mt = mt.filter_cols(hl.literal(set(overlap)).contains(mt.s))
        mt2 = mt2.filter_cols(hl.literal(set(overlap)).contains(mt2.s))

        # Remove any variants from mt that are already in mt2
        mt = mt.filter_rows(~hl.is_defined(mt2.index_rows(mt.locus, mt.alleles)))
       
        # annotate origin
        mt = mt.annotate_rows(wes = 0)
        mt2 = mt2.annotate_rows(wes = 1)

        # combine the two datasets
        mt = tables.order_cols(mt, mt2)
        mt = mt.union_rows(mt2)

    if checkpoint:
        checkpoint_prefix = out_prefix + "_checkpoint"
        mt = mt.checkpoint(checkpoint_prefix, overwrite = True)
    if extract_samples:
        ht_samples = hl.import_table(extract_samples, no_header=True, key='f0', delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_samples[mt.col_key])) 
    if chrom in "X":
        mt = filter_to_females(mt)
    if ancestry:
        mt = samples.filter_ukb_to_ancestry(mt, ancestry)
    if min_mac:
        mt = filter_min_mac(mt, int(min_mac))
    if missing:
        mt = filter_missing(mt, float(missing))
    if missing or min_mac or chrom in "X":
        mt = io.recalc_info(mt)
    if exclude_trio_parents:
        pids = samples.get_parents_by_fam(mt, ["TRIO"])
        if export_parents:
            mt_parents = mt.filter_cols(hl.literal(pids).contains(mt.s))
            io.export_table(mt_parents, out_prefix + "_parents", out_type)
        mt = mt.filter_cols(~hl.literal(pids).contains(mt.s))
    
    # always export matrix table 
    if out_type != "mt":
        io.export_table(mt, out_prefix, "mt")
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='chromosome')
    parser.add_argument('--input_wes_path', default=None, help='What is the input path to the file?')
    parser.add_argument('--input_calls_path', default=None, help='What is the input path to the file?')
    parser.add_argument('--input_wes_type', default=None, help='What input type?')
    parser.add_argument('--input_calls_type', default=None, help='What input type?')
    parser.add_argument('--exclude_trio_parents', default=None, action='store_true', help='Exclude parents of duo/trio relationships')
    parser.add_argument('--ancestry', default=None, help='filter to specific ancestry')
    parser.add_argument('--convert_sample_id', default=None, action='store_true', help='convert to lindgren sample id')
    parser.add_argument('--export_parents', default=None, action='store_true', help='Export parents genotypes seperately')
    parser.add_argument('--checkpoint', default=None, action='store_true', help='Checkpoint after combining rows')
    parser.add_argument('--extract_samples', default=None, help='HailTable with samples to be extracted.')
    parser.add_argument('--min_mac', default=None, help='Filter to MAC >= value')
    parser.add_argument('--missing', default=None, help='Filter to variants to have le value in genotype missingness')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type out vcf/plink/mt')
    args = parser.parse_args()

    main(args)


