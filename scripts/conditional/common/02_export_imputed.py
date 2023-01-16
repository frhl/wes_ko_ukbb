#!/usr/bin/env python3

import hail as hl
import argparse
import os.path
import sys

from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import variants
from ko_utils import io

def main(args):

    chrom = args.chrom
    extract = args.extract
    extract2 = args.extract2
    min_info = float(args.min_info)
    min_maf = float(args.min_maf)
    missing = float(args.missing)
    out_prefix = args.out_prefix
    out_type = args.out_type

    reference_genome = 'GRCh38'
    hail_init.hail_bmrc_init(log='logs/hail/export_imputed.log', default_reference=reference_genome) 

    checkpoint_dir = out_prefix + "_checkpoint.mt"
    checkpoint_file = checkpoint_dir + "/_SUCCESS"

    if not os.path.isfile(checkpoint_file):
        # load chromosome
        mfi = genotypes.get_ukb_imputed_v3_mfi(chrom)
        mfi = mfi.annotate(chrom = chrom)

        # annotate mfi
        mfi = mfi.annotate(ref = hl.if_else(mfi.f6 == mfi.a1, mfi.a2, mfi.a1))
        mfi = mfi.annotate(variant = hl.delimit([hl.str(mfi.chrom), hl.str(mfi.position), mfi.ref, mfi.f6], ':'))
        mfi = mfi.key_by(**hl.parse_variant(mfi.variant,  reference_genome= 'GRCh37'))

        # load imputed
        mt = genotypes.get_ukb_imputed_v3_bgen([chrom])
        mt = mt.repartition(64)
        mt = mt.annotate_rows(info_score = mfi[mt.row_key].info)
        mt = mt.select_entries(mt.GT)
        
        # Filter to relevant samples
        if extract:
            ht_final_samples = hl.import_table(
                extract, no_header=True, key='f0', delimiter=',')
            mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))

        # Filter to relevant samples
        if extract2:
            ht_final_samples2 = hl.import_table(
                extract2, no_header=True, key='f0', delimiter=',')
            mt = mt.filter_cols(hl.is_defined(ht_final_samples2[mt.col_key]))

        # perform variant filtering after subsetting samples
        info_expr = mt.info_score > min_info
        maf_expr = variants.get_maf_expr(mt) > min_maf
        missing_expr = variants.get_missing_expr(mt) < missing
        mt = mt.filter_rows((info_expr) & (maf_expr) & (missing_expr))
        mt = mt.repartition(64)
        mt = mt.checkpoint(out_prefix + "_checkpoint.mt", overwrite = True)
    else:
        mt = io.import_table(checkpoint_dir, "mt", calc_info = False)
        mt = mt.repartition(64)

    # perform liftover
    mt = variants.liftover(mt, fix_ref = False)
    
    # add variant IDs
    mt = mt.annotate_rows(
                varid = hl.delimit(
                    [hl.str(mt.locus.contig),
                     hl.str(mt.locus.position),
                     mt.alleles[0],
                     mt.alleles[1]],
                    ':')
                ) 

    # export imputed
    io.export_table(mt, out_prefix, out_type)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, required=True, help='chromosome to load')
    parser.add_argument('--extract', default=None, help='Path to HailTable that contains the final samples included in the analysis.')
    parser.add_argument('--extract2', default=None, help='Path to HailTable that contains the final samples included in the analysis.')
    parser.add_argument('--min_maf', default=0.01, help='What min_maf threshold should be used?')
    parser.add_argument('--min_info', default=0.8, help='What info threshold should be used?')
    parser.add_argument('--missing', default=0.1, help='What info threshold should be used?')
    parser.add_argument('--out_prefix', default=None, required=True, help='Path prefix for output dataset (plink format)')
    parser.add_argument('--out_type', default=None, required=True, help='Path prefix for output dataset (plink format)')
    args = parser.parse_args()

    main(args)


