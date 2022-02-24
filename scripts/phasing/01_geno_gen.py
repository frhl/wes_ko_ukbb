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

from ko_utils.variants import filter_min_mac, filter_missing

def main(args):

    # parser
    chrom = args.chrom
    input_path = args.input_path
    input_type = args.input_type
    chrom = args.chrom
    dataset = args.dataset
    convert_sample_id = args.convert_sample_id
    min_info = args.min_info
    liftover = args.liftover
    min_mac = args.min_mac
    missing = args.missing
    ancestry = args.ancestry
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    
    if dataset in "imp":
        mt = genotypes.get_ukb_imputed_v3_bgen(chroms=[chrom])
        if min_info:
            ht = genotypes.get_ukb_parsed_imputed_v3_mfi(chroms=[chrom])
            mt = mt.annotate_rows(info_score = ht[(mt.locus, mt.alleles)].info)
            mt = mt.filter_rows(mt.info_score >= float(min_info))
    elif dataset in "calls":
        mt = genotypes.get_ukb_genotypes_bed(chroms=[chrom])
    else:
        raise TypeError(f"{dataset} is not 'imp' or 'calls'")

    if convert_sample_id:
        mt = samples.convert_sample_ids(mt, 12788, 11867)
    if liftover:
        mt = variants.liftover(mt, from_build='GRCh37', to_build='GRCh38', drop_annotations=True)

    if input_path and input_type:

        mt2 = qc.get_table(input_path, input_type) # assuming build GRCh38
        mt2_count = mt2.count()
        mt_count = mt.count()
        print(f"mt2 count: {mt2_count} .. mt1 count: {mt_count}")
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
        mt2_count = mt2.count()
        mt_count = mt.count()
        print(f"again .. mt2 count: {mt2_count} .. mt1 count: {mt_count}")
        
        # combine the two datasets
        mt = tables.order_cols(mt, mt2)
        mt = mt.union_rows(mt2)
        final = mt.count()
        print(f"Final aggregated count: {final}")
    
    if missing or min_mac:
        mt = io.recalc_info(mt)
    if ancestry:
        mt = samples.filter_ukb_to_ancestry(mt, ancestry)
    if min_mac:
        mt = filter_min_mac(mt, int(min_mac))
    if missing:
        mt = filter_missing(mt, float(missing))
   
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='chromosome')
    parser.add_argument('--input_path', default=None, help='What is the input path to the file?')
    parser.add_argument('--input_type', default=None, help='What input type?')
    parser.add_argument('--min_info', default=None, help='minimum info score (only imputed data)')
    parser.add_argument('--liftover', default=None, action='store_true', help='perform liftover')
    parser.add_argument('--ancestry', default=None, help='filter to specific ancestry')
    parser.add_argument('--convert_sample_id', default=None, action='store_true', help='convert to lindgren sample id')
    parser.add_argument('--dataset', default=None, help='Either "imp" or "calls".')
    parser.add_argument('--min_mac', default=None, help='Filter to MAC >= value')
    parser.add_argument('--missing', default=None, help='Filter to variants to have le value in genotype missingness')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type out vcf/plink/mt')
    args = parser.parse_args()

    main(args)


