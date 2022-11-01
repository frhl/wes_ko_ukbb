#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import variants
from ukb_utils import samples
from ko_utils.samples import filter_to_females
from ko_utils.variants import filter_min_mac, filter_missing

def main(args):

    # parser
    chrom = args.chrom
    dataset = args.dataset
    convert_sample_id = args.convert_sample_id
    extract_samples = args.extract_samples
    filter_incorrect_reference = args.filter_incorrect_reference
    exclude_trio_parents = args.exclude_trio_parents
    export_parents = args.export_parents
    min_info = args.min_info
    liftover = args.liftover
    min_mac = args.min_mac
    min_maf = args.min_maf
    missing = args.missing
    ancestry = args.ancestry
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/01_geno_gen.log', 'GRCh38')

    if dataset in "imp":
        mt = genotypes.get_ukb_imputed_v3_bgen(chroms=[chrom])
        if min_info:
            ht = genotypes.get_ukb_parsed_imputed_v3_mfi(chroms=[chrom])
            mt = mt.annotate_rows(info_score = ht[(mt.locus, mt.alleles)].info)
            mt = mt.filter_rows(mt.info_score >= float(min_info))
    elif dataset in "calls":
        mt = genotypes.get_ukb_genotypes_bed(chroms=[chrom])
    else:
        raise TypeError(f"{dataset} is not 'imp' or 'calls'")

    if convert_sample_id:
        mt = samples.convert_sample_ids(mt, 12788, 11867)
    if extract_samples:
        ht_samples = hl.import_table(extract_samples, no_header=True, key='f0', delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_samples[mt.col_key])) 
    if chrom in "X":
        mt = filter_to_females(mt)
    if liftover:
        mt = variants.liftover(mt, from_build='GRCh37', to_build='GRCh38', drop_annotations=True, fix_ref=True)
        mt = mt.filter_rows(mt.locus.contig == "chr" + str(chrom))
    if filter_incorrect_reference:
        mismatch_expr = variants.get_reference_mismatch_expr(mt.locus, mt.alleles, "GRCh38")
        mismatch_n = mt.aggregate_rows(hl.agg.sum(mismatch_expr))
        if mismatch_n > 0:
            print(f"Found {mismatch_n} variants with reference allele not in fasta!")
            mt = mt.filter_rows(~mismatch_expr)
    if ancestry:
        mt = samples.filter_ukb_to_ancestry(mt, ancestry)
    if min_mac:
        mt = mt.filter_rows(variants.get_mac_expr(mt) >= int(min_mac))
    if min_maf:
        mt = mt.filter_rows(variants.get_maf_expr(mt) >= float(min_maf))
    if missing:
        mt = filter_missing(mt, float(missing))
    if exclude_trio_parents:
        mt = mt.checkpoint(out_prefix + ".mt")
        pids = samples.get_parents_by_fam(mt, ["TRIO"])
        if export_parents:
            mt_parents = mt.filter_cols(hl.literal(pids).contains(mt.s))
            io.export_table(mt_parents, out_prefix + "_parents", out_type)
        mt = mt.filter_cols(~hl.literal(pids).contains(mt.s))
    mt = io.recalc_info(mt)

    # always export matrix table 
    if exclude_trio_parents:
        io.export_table(mt, out_prefix + "_no_parents", out_type)
    else:
        io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='chromosome')
    parser.add_argument('--min_info', default=None, help='minimum info score (only imputed data)')
    parser.add_argument('--liftover', default=None, action='store_true', help='perform liftover')
    parser.add_argument('--ancestry', default=None, help='filter to specific ancestry')
    parser.add_argument('--convert_sample_id', default=None, action='store_true', help='convert to lindgren sample id')
    parser.add_argument('--filter_incorrect_reference', default=None, action='store_true', help='Remove any sites with incorrect reference.')
    parser.add_argument('--exclude_trio_parents', default=None, action='store_true', help='Exclude parents of duo/trio relationships')
    parser.add_argument('--export_parents', default=None, action='store_true', help='Export parents genotypes seperately')
    parser.add_argument('--dataset', default=None, help='Either "imp" or "calls".')
    parser.add_argument('--extract_samples', default=None, help='HailTable with samples to be extracted.')
    parser.add_argument('--min_mac', default=None, help='Filter to MAC >= value')
    parser.add_argument('--min_maf', default=None, help='Filter to MAF >= value')
    parser.add_argument('--missing', default=None, help='Filter to variants to have le value in genotype missingness')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type out vcf/plink/mt')
    args = parser.parse_args()

    main(args)


