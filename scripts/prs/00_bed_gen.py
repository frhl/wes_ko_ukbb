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
    exclude_samples = args.exclude_samples
    write_samples = args.write_samples
    min_info = args.min_info
    liftover = args.liftover
    hapmap = args.hapmap
    ancestry = args.ancestry
    dbsnp = args.dbsnp
    filter_missing = args.filter_missing
    random_samples = args.random_samples
    random_seed = args.random_seed
    min_maf = args.min_maf
    only_valid_contigs = args.only_valid_contigs
    out_prefix = args.out_prefix
    out_type = args.out_type
    filter_to_unrelated_using_kinship_coef = args.filter_to_unrelated_using_kinship_coef

    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
   
    if chrom in "AUTOSOMES":
        chrom = list(map(str, range(1, 23)))
    else:
        chrom = [chrom]
        
    if dataset in "imp":
        mt = genotypes.get_ukb_imputed_v3_bgen(chroms=chrom)
        if min_info:
            ht = genotypes.get_ukb_parsed_imputed_v3_mfi(chroms=chrom)
            mt = mt.annotate_rows(info_score = ht[(mt.locus, mt.alleles)].info)
            mt = mt.filter_rows(mt.info_score >= float(min_info))
    elif dataset in "calls":
        mt = genotypes.get_ukb_genotypes_bed(chroms=chrom)
    else:
        raise TypeError(f"{dataset} is not 'imp' or 'calls'")
    
    if extract_samples:
        ht_samples = hl.import_table(extract_samples, no_header=True, key='f0', delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_samples[mt.col_key]))
    
    if exclude_samples:
        ht_samples = hl.import_table(exclude_samples, no_header=True, key='f0', delimiter=',')
        mt = mt.filter_cols(~hl.is_defined(ht_samples[mt.col_key]))
  
    if filter_to_unrelated_using_kinship_coef:
        mt = samples.filter_ukb_to_unrelated_using_kinship(mt)

    if ancestry:
        mt = samples.filter_ukb_to_ancestry(mt, ancestry)

    if random_samples:
        mt = samples.choose_col_subset(mt, int(random_samples), seed = int(random_seed))

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

    if dbsnp:
        ht = variants.get_dbsnp_table(version=155, build='GRCh38')
        mt = mt.annotate_rows(rsid = ht.rows()[mt.row_key].rsid)

    if only_valid_contigs:
        chroms = [f'chr{x}' for x in range(1,23)]
        chroms.append('chrX')
        filter_expr = hl.literal(set(chroms)).contains(mt.locus.contig)
        mt = mt.filter_rows(filter_expr)
    
    if write_samples and out_prefix:
        mt.cols().write(out_prefix + "_samples.ht", overwrite = True)

    if out_type and out_prefix:
        io.export_table(mt, out_prefix, out_type)


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='chromosome')
    parser.add_argument('--min_info', default=None, help='minimum info score (only imputed data)')
    parser.add_argument('--liftover', default=None, action='store_true', help='perform liftover')
    parser.add_argument('--dbsnp', default=None, action='store_true', help='Annotate rsids.')
    parser.add_argument('--filter_to_unrelated_using_kinship_coef', default=None, action='store_true', help='Exclude any related individuals.')
    parser.add_argument('--filter_missing', default=None, help='Filter to variants with lt value in genotype missingness.')
    parser.add_argument('--extract_samples', default=None, help='Subset to sample IDs in MatrixTable')
    parser.add_argument('--exclude_samples', default=None, help='Exclude sample IDs from MatrixTable')
    parser.add_argument('--write_samples', default=None, action='store_true', help='Write hail table cols of subsetted dataset')
    parser.add_argument('--dataset', default=None, help='Either "imp" or "calls".')
    parser.add_argument('--ancestry', default=None, help='Either "eur" or "all".')
    parser.add_argument('--hapmap', default=None, help='Path to hapmap file')
    parser.add_argument('--random_samples', default=None, help='Subset to random samples')
    parser.add_argument('--random_seed', default=42, help='Seed for randomizer')
    parser.add_argument('--only_valid_contigs', default=None, action='store_true', help='Subset variants only normal contigs chr1..22x')
    parser.add_argument('--min_maf', default=None, help='Subset to variants based on minimum MAF')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Either "mt", "vcf" or "plink"')
    args = parser.parse_args()

    main(args)


