#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import variants
from ukb_utils import samples


def main(args):

    # parser
    chrom = args.chrom
    dataset = args.dataset
    extract_samples = args.extract_samples
    min_info = args.min_info
    liftover = args.liftover
    hapmap = args.hapmap
    ancestry = args.ancestry
    dbsnp = args.dbsnp
    filter_missing = args.filter_missing
    exclude_related = args.exclude_related
    random_samples = args.random_samples
    random_variants = args.random_variants
    min_maf = args.min_maf
    input_annotation_path = args.input_annotation_path
    sample_seed = args.sample_seed
    variant_seed = args.variant_seed
    in_prefix = args.in_prefix
    in_type = args.in_type
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
   
    if chrom in "AUTOSOMES":
        chrom = list(map(str, range(1, 23)))
    else:
        chrom = [chrom]
    
    if dataset:
        if dataset in "imp":
            mt = genotypes.get_ukb_imputed_v3_bgen(chroms=chrom)
            if min_info:
                ht = genotypes.get_ukb_parsed_imputed_v3_mfi(chroms=chrom)
                mt = mt.annotate_rows(info_score = ht[(mt.locus, mt.alleles)].info)
                mt = mt.filter_rows(mt.info_score >= float(min_info))
        elif dataset in "calls":
            mt = genotypes.get_ukb_genotypes_bed(chroms=chrom)
    else:
        mt = io.import_table(in_prefix, in_type)
        print(mt.count())

    if extract_samples:
        ht_samples = hl.import_table(extract_samples, no_header=True, key='f0', delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_samples[mt.col_key]))
  
    if exclude_related:
        related = samples.get_ukb_is_related_using_kinship_expr(mt)
        mt = mt.filter_cols(~related) 

    if ancestry:
        mt = samples.filter_ukb_to_ancestry(mt, ancestry)

    if random_samples:
        mt = samples.choose_col_subset(mt, int(random_samples), seed = int(sample_seed))

    #if random_variants:
    #    mt = variants.choose_row_subset(mt, int(random_variants), seed = int(variant_seed))

    if hapmap:
        ht = hl.read_table(hapmap)
        ht = ht.key_by('locus_grch37')
        mt = mt.filter_rows(hl.is_defined(ht[mt.locus]))
    
    if filter_missing:
        missing = hl.agg.mean(hl.is_missing(mt.GT)) <= float(filter_missing)
        mt = mt.filter_rows(missing)

    if min_maf:
        mt = mt.filter_rows(variants.get_maf_expr(mt) > float(min_maf))

    if liftover:
        mt = variants.liftover(mt, from_build='GRCh37', to_build='GRCh38', drop_annotations=True)
    
    if input_annotation_path:
        consequence_annotations = hl.read_table(input_annotation_path)
        mt = mt.annotate_rows(consequence=consequence_annotations[mt.row_key]) 
    
    if dbsnp:
        ht = variants.get_dbsnp_table(version=155, build='GRCh38')
        mt = mt.annotate_rows(rsid = ht.rows()[mt.row_key].rsid)

    io.export_table(mt, out_prefix, out_type)


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='chromosome')
    parser.add_argument('--min_info', default=None, help='minimum info score (only imputed data)')
    parser.add_argument('--liftover', default=None, action='store_true', help='perform liftover')
    parser.add_argument('--dbsnp', default=None, action='store_true', help='Annotate rsids.')
    parser.add_argument('--input_annotation_path', default=None, help='Annotate variant consequence using a hail table')
    parser.add_argument('--exclude_related', default=None, action='store_true', help='Exclude any related individuals.')
    parser.add_argument('--filter_missing', default=None, help='Filter to variants with lt value in genotype missingness.')
    parser.add_argument('--extract_samples', default=None, help='Subset to sample IDs in file')
    parser.add_argument('--dataset', default=None, help='Either "imp" or "calls".')
    parser.add_argument('--ancestry', default=None, help='Either "eur" or "all".')
    parser.add_argument('--hapmap', default=None, help='Path to hapmap file')
    parser.add_argument('--random_samples', default=None, help='Subset to random samples')
    parser.add_argument('--random_variants', default=None, help='Subset to random variants')
    parser.add_argument('--sample_seed', default=None, help='What seed should be used for randomization?')
    parser.add_argument('--variant_seed', default=None, help='What seed should be used for randomization?')
    parser.add_argument('--min_maf', default=None, help='Subset to variants based on minimum MAF')
    parser.add_argument('--in_prefix', default=None, help='Path prefix for input dataset')
    parser.add_argument('--in_type', default=None, help='Either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Either "mt", "vcf" or "plink"')
    args = parser.parse_args()

    main(args)


