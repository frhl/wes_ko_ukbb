#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import variants


def main(args):

    filter_invariant = args.filter_invariant
    common_path = args.common_path
    common_type = args.common_type
    write_samples = args.write_samples
    hapmap = args.hapmap
    dbsnp = args.dbsnp
    filter_missing = args.filter_missing
    min_maf = args.min_maf
    only_valid_contigs = args.only_valid_contigs
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
   
    mt = io.import_table(common_path, common_type)
    
    if hapmap:
        ht = hl.read_table(hapmap)
        ht = ht.key_by(*[ht.locus, ht.alleles])
        ht = ht.key_by(**hl.parse_variant(ht.grch38_varid, reference_genome = 'GRCh38'))
        mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))
    
    if filter_missing:
        missing = hl.agg.mean(hl.is_missing(mt.GT)) <= float(filter_missing)
        mt = mt.filter_rows(missing)

    if min_maf:
        mt = mt.filter_rows(variants.get_maf_expr(mt) > float(min_maf))

    if filter_invariant:
        mt = mt.annotate_rows(stdev = hl.agg.stats(mt.GT.n_alt_alleles()).stdev)
        expr_invariant = mt.stdev > 0
        mt = mt.filter_rows(expr_invariant)

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

    # required for match
    mt = mt.transmute_rows(
            rsid = hl.delimit(
                [hl.str(mt.locus.contig), 
                 hl.str(mt.locus.position), 
                 hl.str(mt.alleles[0]), 
                 hl.str(mt.alleles[1])], ':')
            )    

    if out_type and out_prefix:
        io.export_table(mt, out_prefix, out_type)


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--dbsnp', default=None, action='store_true', help='Annotate rsids.')
    parser.add_argument('--filter_missing', default=None, help='Filter to variants with lt value in genotype missingness.')
    parser.add_argument('--common_path', default=None, help='path to processed common variants that have been imputed.')
    parser.add_argument('--common_type', default="mt", help='path to processed common variants that have been imputed.')
    parser.add_argument('--write_samples', default=None, action='store_true', help='Write hail table cols of subsetted dataset')
    parser.add_argument('--hapmap', default=None, help='Path to hapmap file')
    parser.add_argument('--only_valid_contigs', default=None, action='store_true', help='Subset variants only normal contigs chr1..22x')
    parser.add_argument('--filter_invariant', default=None, action='store_true', help='Subset variants only normal contigs chr1..22x')
    parser.add_argument('--min_maf', default=None, help='Subset to variants based on minimum MAF')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Either "mt", "vcf" or "plink"')
    args = parser.parse_args()

    main(args)


