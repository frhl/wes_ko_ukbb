#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import qc

def main(args):
    
    # parser
    input_path = args.input_path
    input_type = args.input_type
    initial_variant_list = args.initial_variant_list

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')

    
    ukb_bed_X="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chrX_v2.bed"
    ukb_bim_X="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chrX_v2.bim"
    ukb_bed_Y="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chrY_v2.bed"
    ukb_bim_Y="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chrY_v2.bim"
    ukb_fam="/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"
    INITIAL_SAMPLES = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/03_initial_qc.keep.sample_list'
    PHENOTYPES_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/' + TRANCHE + '/phenotypes.ht'
    PRUNED_CHRX_VARIANTS = '/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_chrX_pruned.prune.in'

    # Outputs
    IMPUTESEX_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/04_imputesex.ht'
    IMPUTESEX_FILE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/04_imputesex.tsv.bgz'
    Y_NCALLED = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/04_ycalled.tsv.bgz'

    ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
    ht_pruned_chrx_variants = hl.import_table(PRUNED_CHRX_VARIANTS, no_header=True)
    ht_pruned_chrx_variants = ht_pruned_chrx_variants.transmute(rsid=ht_pruned_chrx_variants.f0).key_by('rsid')
    sample_annotations = hl.read_table(PHENOTYPES_TABLE)

    # Read in the plink file
    mt = hl.import_plink(bed=ukb_bed_X, bim=ukb_bim_X, fam=ukb_fam, reference_genome='GRCh37')
    mt = mt.key_rows_by(mt.rsid)
    mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
    mt = mt.filter_rows(hl.is_defined(ht_pruned_chrx_variants[mt.row_key]))
    mt = mt.key_rows_by(mt.locus, mt.alleles)
    n = mt.count()

    print('n samples:')
    print(n[1])
    print('n variants:')
    print(n[0])

    imputed_sex = hl.impute_sex(mt.GT, female_threshold=0.6, male_threshold=0.6)
    mt = mt.annotate_cols(phenotype = sample_annotations[mt.s])
    mt = mt.annotate_cols(impute_sex = imputed_sex[mt.s])

    mt.cols().select('impute_sex', 'phenotype').flatten().export(IMPUTESEX_FILE)
    mt.cols().write(IMPUTESEX_TABLE, overwrite=True)

    # Now, look on the Y chromosome, and determine non-missing allele count on the y.
    mt = hl.import_plink(bed=ukb_bed_Y, bim=ukb_bim_Y, fam=ukb_fam, reference_genome='GRCh37')
    mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
    mt = mt.filter_rows(mt.locus.in_y_nonpar() | mt.locus.in_y_par())
    mt = hl.sample_qc(mt, name='qc')

    mt_cols = mt.cols()
    mt_cols.select(n_called=mt_cols.qc.n_called).export(Y_NCALLED)



if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--initial_variant_list', default=None, help='List of initial QCed variants')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
 
    args = parser.parse_args()

    main(args)

