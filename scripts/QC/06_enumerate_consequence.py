#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import qc
from ko_utils import analysis
from gnomad.utils.vep import process_consequences

def main(args):
    
    # parser
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
   
    vep_path_list = args.vep_path_list
    initial_samples_list = args.initial_samples_list
    initial_variants_list = args.initial_variants_list
    sexcheck_list = args.sexcheck_list

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    #mt = qc.get_table(input_path=input_path, input_type=input_type)

    # Inputs
    #MT_HARDCALLS = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_hardcalls_chr' + CHR + '.mt'
    #UKB_vep_output = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_variants_vep/' + TRANCHE + '/'
    #ANNOTATION_TABLE = UKB_vep_output + 'ukb_wes_' + TRANCHE + '_filtered_chr' + CHR + '_vep_qc.ht'
    #PHENOTYPES_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/' + TRANCHE + '/phenotypes.ht'
    #INITIAL_SAMPLES = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/03_initial_qc.keep.sample_list'
    #INITIAL_VARIANT_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/variants/02_prefilter_chr' + CHR +'.keep.variant.ht'
    #SEXCHECK_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/04_sexcheck.remove.sample_list'

    # Outputs
    #URV_FILE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/07_URVs_chr' + CHR + '.tsv.bgz'

    ht_initial_variants = hl.read_table(initial_variants_list)
    ht_initial_samples = hl.import_table(initial_samples_list, no_header=True, key='f0')
    ht_sexcheck_samples = hl.import_table(sexcheck_list, no_header=True, key='f0')

    # Annotate variants with counts from non-psychiatric version of Gnomad.
    # Fill in variants not in Gnomad variant list.
    # Annotate variants with LoF/damaging missense annotation.

    mt = qc.get_table(input_path, input_type=input_type)
    mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
    mt = mt.filter_cols(~hl.is_defined(ht_sexcheck_samples[mt.col_key]))
    mt = mt.filter_rows(hl.is_defined(ht_initial_variants[mt.row_key]))

    # Drop some fields that are not needed.
    mt = mt.drop('rsid', 'info', 'filters')

    #sample_annotations = hl.read_table(PHENOTYPES_TABLE)
    consequence_annotations = hl.read_table(ANNOTATION_TABLE)
    #output_vep_ht_path = UKB_vep_output + 'ukb_wes_' + TRANCHE + '_filtered_chr' + CHR + '_vep.ht'

    #mt = mt.annotate_cols(phenotype = sample_annotations[mt.s])
    mt = mt.annotate_rows(consequence = consequence_annotations[mt.row_key])
    mt = mt.annotate_rows(is_singleton = hl.agg.sum(mt.GT.n_alt_alleles()) == 1)
    #mt = mt.filter_rows((mt.is_singleton) & (~mt.consequence.inGnomAD))

    mt = mt.annotate_cols(n_coding_URV_SNP = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_snp(mt.alleles[0], mt.alleles[1]) & (mt.vep.consequence_category != "non_coding")),
                      n_coding_URV_indel = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_indel(mt.alleles[0], mt.alleles[1]) & (mt.vep.consequence_category != "non_coding")),
                      n_URV_PTV = hl.agg.count_where(mt.GT.is_non_ref() & (mt.consequence_category == "ptv")),
                      n_URV_damaging_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.vep.consequence_category == "damaging_missense")),
                      n_URV_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.vep.consequence_category == "missense")),
                      n_URV_synonymous = hl.agg.count_where(mt.GT.is_non_ref() & (mt.vep.consequence_category == "synonymous")))

    mt.cols().flatten().export(out_prefix + '.tsv.bgz')


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--vep_path_list', default=None, help='Path to VEP file with dbNSFP hail table"')
    parser.add_argument('--initial_samples_list', default=None, help='Path to initial QCed samples hail table')
    parser.add_argument('--initial_variants_list', default=None, help='Path to initial QCed variant hail table')
    parser.add_argument('--sexcheck_list', default=None, help='Path to sexcheck (imputed sex)')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')


    args = parser.parse_args()

    main(args)

