#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import samples
from ko_utils import qc
from ko_utils import analysis


def main(args):

    # parser
    chrom = args.chrom
    input_annotation_path = args.input_annotation_path
    out_prefix = args.out_prefix

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')  # from zulip

    # load qc table and annotate
    hts = list()
    for CHR in range(1,22):

        # get paths
        QC_PATH  = "data/qc_old/ukb_wes_200k_chr" + str(CHR) + "_variants_unphased.ht"
        ANNOTATION_TABLE = "data/vep/hail/ukb_wes_200k_chr" + str(CHR) + "_vep.ht"
        
        # load files
        ht = hl.read_table(QC_PATH)
        vep = hl.read_table(ANNOTATION_TABLE)
        
        # annotate table
        ht = ht.annotate(worst_csq_for_variant_canonical=vep[ht.key].vep.worst_csq_for_variant_canonical)
        annotations = analysis.annotation_case_builder(ht.worst_csq_for_variant_canonical, use_loftee=True)
        ht = ht.annotate(consequence_category=annotations)

        hts += ht

    # combine all rows
    ht = hts[0].union(*hts[1:])
    
    # combine by category
    ht_cat_all = ht.group_by(ht.consequence_category).aggregate(n=hl.agg.count())
    ht_cat_all.flatten().export(out_prefix + '_category.tsv')

    # get AC per individual
    quantile = [0, 0.25,0.50,0.75, 1]
    ht_cat_quantile = ht.group_by(ht.consequence_category).aggregate(n=hl.agg.approx_quantiles(ht.variant_qc.AC[1], quantiles))
    ht_cat_quantile.flatten().export(out_prefix + '_af_quantiles.tsv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_qc_path', default=None, help='Path to input')
    parser.add_argument('--input_annotation_path',default=None,help='path to HailTable with VEP and dbNSFP annotations')
    parser.add_argument('--out_prefix',default=None,help='Path prefix for output dataset')
    parser.add_argument('--chrom', default=None, help='Chromosome to be used')

    args = parser.parse_args()

    main(args)
