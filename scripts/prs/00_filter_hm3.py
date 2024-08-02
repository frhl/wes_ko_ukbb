#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import variants
from ukb_utils import genotypes


def main(args):

    chrom = args.chrom
    hapmap = args.hapmap
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip

    # get chromosomes 
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

    # subset to hapmap3 on reference GRCH37
    ht = hl.read_table(hapmap)
    ht = ht.key_by(**hl.parse_variant(ht.grch37_varid, reference_genome = 'GRCh37'))
    mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))

    # liftover remaining variants
    mt = variants.liftover(mt, fix_ref = False)
    mt = mt.repartition(16)

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

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='Path to hapmap file')
    parser.add_argument('--hapmap', default=None, help='Path to hapmap file')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Either "mt", "vcf" or "plink"')
    args = parser.parse_args()

    main(args)


