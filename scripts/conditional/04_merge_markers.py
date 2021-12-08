import os
import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import samples
from ukb_utils import variants
from ko_utils import qc


def main(args):
   
    input_vcf = args.input_vcf
    input_markers = args.input_markers
    out_prefix = args.out_prefix

    reference_genome = 'GRCh38'
    hail_init.hail_bmrc_init_local(
        'logs/hail/hail_conditional_tables.log',
        reference_genome)
    hl._set_flags(no_whole_stage_codegen='1')

    # load table of genes that are significant in primary analysis
    ht = hl.import_table(input_markers)

    # table of variants to condition on
    ht = ht.annotate(variant = hl.delimit([ht.f0, ht.f1, ht.f3, ht.f4], ':'))
    ht = ht.key_by(**hl.parse_variant(ht.variant, reference_genome = 'GRCh37'))

    # get genotype for selected variants
    chromosomes = [x.replace('chr', '') for x in list(set(ht.locus.contig.collect()))]
    mtc = genotypes.get_ukb_imputed_v3_bgen(chromosomes)
    mtc = variants.liftover(mtc)
    mtc = mtc.filter_rows(hl.is_defined(ht[mtc.row_key]))
    mtc = mtc.annotate_entries(DS = mtc.GT.n_alt_alleles())

    # load knockout file
    mt = hl.import_vcf(input_vcf)
    
    # order mtc and mt

    
    
    
    #hl.export_vcf(mt, out_prefix + '.vcf.bgz')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_vcf', default=None, help='path to VCF file')
    parser.add_argument('--input_markers', default=None, help='Table of markers to be inputted')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset (plink format)')
    args = parser.parse_args()

    main(args)
