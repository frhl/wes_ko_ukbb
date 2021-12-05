import os
import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import samples
from ukb_utils import variants
from ko_utils import qc


def main(args):

    out_prefix = args.out_prefix
    gene_table = args.gene_table
    final_sample_list = args.final_sample_list
    phenotype = args.phenotype
    annotation = args.annotation
    padding = int(args.padding)

    reference_genome = 'GRCh38'
    hail_init.hail_bmrc_init_local(
        'logs/hail/hail_conditional_tables.log',
        reference_genome)
    hl._set_flags(no_whole_stage_codegen='1')

    # load table of genes that are significant in primary analysis

    # Get chromosomes
    chromosomes = [x.replace('chr', '') for x in list(set(ht.contig.collect()))]
    mt = genotypes.get_ukb_genotypes_bed(chromosomes)
    mt = mt.annotate_rows(**{'info': hl.agg.call_stats(mt.GT, mt.alleles)})

    # filter samples
    if final_sample_list:
        ht_final_samples = hl.import_table(
            final_sample_list, no_header=True, key='f0', delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))

    # filter variants based on missingness
    mt = qc.filter_min_missing(mt, 0.10)
    mt = mt.annotate_rows(info=mt.info.annotate(AC=mt.info.AC[1]))
    mt = mt.annotate_rows(info=mt.info.annotate(AF=mt.info.AF[1]))
    mt = qc.filter_min_maf(mt, 0.01)

    # perform liftover to GRCh38 and filter to intervals
    mt = variants.liftover(mt)
    hail_intervals = ht.intervals.collect()
    mt = hl.filter_intervals(mt, hail_intervals)
    print(hail_intervals)

    n = mt.count()
    print(f'filtered genotypes to {n} variants/samples. Writing to VCF..')
    hl.export_vcf(mt, out_prefix + '.vcf.bgz')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--padding', default=0, help='How much extra padding should be included around genes (upstream/downsream)?')
    parser.add_argument('--final_sample_list', default=None, help='Path to HailTable that contains the final samples included in the analysis.')
    parser.add_argument('--gene_table', default=None, help='Path to HailTable that contains the significant genes from primary analysis.')
    parser.add_argument('--phenotype', default=None, help='Path to phenotype table that also includes sex.')
    parser.add_argument('--annotation', default=None, help='String. What mutation subset should be applied?')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset (plink format)')
    args = parser.parse_args()

    main(args)
