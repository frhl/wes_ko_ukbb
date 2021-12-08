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

    # load knockout VCF
    mt = hl.import_vcf(input_vcf)

    # Subset to conditioning markers
    geno = genotypes.get_ukb_genotypes_bed(AUTOSOMES)
    geno = geno.filter_cols(hl.is_defined(geno[mt.row_key].s))
    geno = variants.liftover(geno)
    geno = geno.filter_rows(hl.is_defined(geno[ht.key]))
    n = geno.count()
    print(f'Filtered to {n} genotyped variants. Merging..')
    
    #mt.union_rows(geno)
    
    hl.export_vcf(mt, out_prefix + '.vcf.bgz')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_vcf', default=None, help='path to VCF file')
    parser.add_argument('--input_markers', default=None, help='Table of markers to be inputted')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset (plink format)')
    args = parser.parse_args()

    main(args)
