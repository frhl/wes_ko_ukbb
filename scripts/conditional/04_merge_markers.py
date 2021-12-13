import os
import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import variants
from ukb_utils import tables


def main(args):
   
    input_mt = args.input_mt
    input_markers = args.input_markers
    chrom = args.chrom
    out_prefix = args.out_prefix

    reference_genome = 'GRCh38'
    hail_init.hail_bmrc_init_local(
        'logs/hail/hail_conditional_tables.log',
        reference_genome)
    hl._set_flags(no_whole_stage_codegen='1')

    # load table of genes that are significant in primary analysis
    ht = hl.import_table(input_markers, no_header = True, delimiter = ' ')

    # table of variants to condition on
    ht = ht.annotate(variant = hl.delimit([ht.f0, ht.f1, ht.f3, ht.f4], ':'))
    ht = ht.key_by(**hl.parse_variant(ht.variant, reference_genome = 'GRCh38'))
    ht = ht.filter(hl.literal(chrom).contains(ht.locus.contig))
    rows = int(ht.count())

    if rows > 0:

        # get genotype for selected variants
        chromosomes = [x.replace('chr', '') for x in list(set(ht.locus.contig.collect()))]
        imp = genotypes.get_ukb_imputed_v3_bgen(chromosomes)
        imp = variants.liftover(imp)
        imp = imp.filter_rows(hl.is_defined(ht[imp.row_key]))
        imp = imp.annotate_entries(DS = imp.dosage)
        imp = imp.drop('varid','new_locus','new_alleles','old_locus')
        imp = imp.select_entries(imp.DS)

        # load knockout file
        wes = hl.read_matrix_table(input_mt)
        wes = wes.drop('eur')
        n_before = wes.count()
        print(f"Pre-merging variants/samples {n_before}")

        # Get intersecting individuals
        wes_sids = wes.s.collect()
        imp_sids = imp.s.collect()
        overlap = list(set(wes_sids) & set(imp_sids))[1:5]
        wes = wes.filter_cols(hl.literal(set(overlap)).contains(wes.s))
        imp = imp.filter_cols(hl.literal(set(overlap)).contains(imp.s))
        
        # combine the variants
        imp = tables.order_cols(imp, wes)
        mt = imp.union_rows(wes)
        n_after = mt.count()
        print(f"Pre-merging variants/samples {n_after}")

        # order imp and mt
        hl.export_vcf(mt, out_prefix + '.vcf.bgz')

        # write file with conditoning markers.
        locus = hl.delimit([imp.locus.contig, hl.str(imp.locus.position)],':')
        alleles = hl.delimit([imp.alleles[0],imp.alleles[1]],'/')
        variant = hl.delimit([locus, alleles],'_')
        saige_variants = ",".join(variant.collect()) 
        
        with open(out_prefix + ".cond_markers", 'a') as outfile:
            outfile.write(saige_variants)
    else:
        print(f"File '${input_markers}' does not contain any markers for chromosome '{chrom}'. Exiting..")

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_mt', default=None, help='path to MatrixTable directory')
    parser.add_argument('--input_markers', default=None, help='Table of markers to be inputted')
    parser.add_argument('--chrom', default=None, help='Chromosome')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset (plink format)')
    args = parser.parse_args()

    main(args)
