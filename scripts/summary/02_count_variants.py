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
    input_qc_path = args.input_qc_path
    input_annotation_path = args.input_annotation_path
    out_prefix = args.out_prefix
    #out_type   = args.out_type

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')  # from zulip

    # load qc table and annotate
    ht = hl.read_table(input_qc_path)
    vep = hl.read_table(input_annotation_path)
    ht = ht.annotate(
        worst_csq_for_variant_canonical=vep[ht.key].vep.worst_csq_for_variant_canonical)

    # annotate with category
    annotations = analysis.annotation_case_builder(
        ht.worst_csq_for_variant_canonical, use_loftee=True)
    ht = ht.annotate(consequence_category=annotations)

    # combine tables
    ht_all = ht.group_by(
        ht.worst_csq_for_variant_canonical.most_severe_consequence).aggregate(
        n=hl.agg.count())
    ht_singletons = ht.filter((ht.variant_qc.AC[1] == 1) | (ht.variant_qc.AC[0] == 1) )
    ht_singletons = ht_singletons.group_by(
        ht_singletons.worst_csq_for_variant_canonical.most_severe_consequence).aggregate(
        n=hl.agg.count())

    ht_all.flatten().export(out_prefix + '.tsv')
    ht_singletons.flatten().export(out_prefix + '_singletons.tsv')

    # annotate by consequence category
    ht_cat_all = ht.group_by(
        ht.consequence_category).aggregate(
        n=hl.agg.count())
    
    ht_cat_singletons = ht.filter((ht.variant_qc.AC[1] == 1) | (ht.variant_qc.AC[0] == 1))
    ht_cat_singletons = ht_cat_singletons.group_by(
        ht_cat_singletons.consequence_category).aggregate(
        n=hl.agg.count())

    ht_cat_all.flatten().export(out_prefix + '_category.tsv')
    ht_cat_singletons.flatten().export(out_prefix + '_category_singletons.tsv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_qc_path', default=None, help='Path to input')
    parser.add_argument('--input_annotation_path',default=None,help='path to HailTable with VEP and dbNSFP annotations')
    parser.add_argument('--out_prefix',default=None,help='Path prefix for output dataset')
    parser.add_argument('--chrom', default=None, help='Chromosome to be used')

    args = parser.parse_args()

    main(args)
