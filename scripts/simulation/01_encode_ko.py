#!/usr/bin/env python3

import hail as hl
import numpy as np
import pandas as pd
import argparse
import random
import string

from ko_utils import io
from ko_utils import ko
from ukb_utils import hail_init
from ukb_utils import variants


class SplitArgs(argparse.Action):
    r"""Method for splitting input csv into a list"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))


def get_tid(length=5):
    r"""method for getting random ID string for alleles"""
    return ''.join(random.choices(string.ascii_uppercase + string.ascii_lowercase,
                                  k=length))


def main(args):

    in_prefix = args.in_prefix
    in_type = args.in_type
    out_prefix = args.out_prefix
    out_type = args.out_type
    chrom = int(args.chrom)

    maf_max = args.maf_max
    maf_min = args.maf_min
    exclude = args.exclude
    use_loftee = args.use_loftee
    csqs_category = args.csqs_category
    exclude_singletons = args.exclude_singletons

    seed=42

    # import table
    hail_init.hail_bmrc_init('logs/hail/simulate_phenotype.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
    hl.set_global_seed(int(seed))
    mt = io.import_table(in_prefix, in_type)
    
    # filter on MAF 
    if maf_max:
        mt = mt.annotate_rows(MAF=variants.get_maf_expr(mt))
        mt = mt.filter_rows(mt.MAF <= float(maf_max))

    if exclude_singletons:
        rows_before = mt.count_rows()
        mt = mt.annotate_rows(MAC=variants.get_mac_expr(mt))
        mt = mt.filter_rows(mt.MAC>0)
        rows_after = mt.count_rows()
        print("singleton exclusion before/after")
        print(rows_before)
        print(rows_after)


    # annotate with variant consequence and collapse by gene
    mt = mt.explode_rows(mt.consequence.vep.worst_csq_by_gene_canonical)
    mt = mt.annotate_rows(
        consequence_category=ko.csqs_case_builder(
                worst_csq_expr=mt.consequence.vep.worst_csq_by_gene_canonical,
                use_loftee=True))

    # filter to variants with SD > 0 and codign variants
    categories = ['pLoF','LC','damaging_missense']
    mt = mt.filter_rows(hl.agg.stats(mt.GT.n_alt_alleles()).stdev>0)
    mt = mt.filter_rows(hl.literal(set(categories)).contains(mt.consequence_category))

    # collapse to gene level
    gene_expr = mt.consequence.vep.worst_csq_by_gene_canonical.gene_id
    mt = ko.aggr_phase_count_by_expr(mt, gene_expr)

    # annotate knockout type
    expr_pko = ko.calc_prob_ko(mt.hom_alt_n, mt.phased, mt.unphased)
    expr_ko = ko.annotate_knockout(mt.hom_alt_n, expr_pko)
    mt = mt.annotate_entries(pKO=expr_pko, knockout=expr_ko)   

    # What haplotypes are affected?
    mt = mt.annotate_entries(
        H1 = 1*((mt.phased.a1 > 0) | (mt.hom_alt_n > 0)),
        H2 = 1*((mt.phased.a2 > 0) | (mt.hom_alt_n > 0)),
    )

    # combine into single dosage (G) matrix 
    mt = mt.annotate_entries(G=hl.int32(mt.H1+mt.H2))
    mt = mt.annotate_rows(**{'stats_add': hl.agg.stats(mt.G)})
    mt = mt.filter_rows(mt.stats_add.stdev > 0)
    mt = mt.annotate_entries(
        G_add_norm=(mt.G-mt.stats_add.mean)/mt.stats_add.stdev
    )
    
    # G_norm but only with alternate alleles
    mt = mt.annotate_entries(
            G_rec = (hl.case().when(mt.G == 2, mt.G).default(0))
    )
    
    # normalize by additive psuedo variant count
    mt = mt.annotate_entries(
            G_rec_norm_by_add = (mt.G_rec - mt.stats_add.mean)/mt.stats_add.stdev
    )
 
    # normalize by recessive pseudo variant count
    mt = mt.annotate_rows(**{'stats_rec': hl.agg.stats(mt.G_rec)})
    mt = mt.annotate_entries(
            G_rec_norm_by_rec = (mt.G_rec - mt.stats_rec.mean)/mt.stats_rec.stdev
    )
    
    # setup sites and alleles
    rows = mt.count()[0]
    locus = [ "chr%s:%s" % (chrom, str(i+1)) for i in range(rows)]
    ref = [ get_tid(4) for i in range(rows)]
    alt = [ get_tid(4) for i in range(rows)]

    # convert to HailTable
    df = pd.DataFrame({'locus':locus, 'ref':ref, 'alt':alt})
    ht = hl.Table.from_pandas(df)
    ht = ht.add_index()
    ht = ht.key_by('idx')

    # annotate knockout matrix
    prob = mt.annotate_entries(DS=mt.pKO * 2)
    prob = prob.annotate_entries(GT=ko.get_gt_from_floor_ds(prob.DS))
    prob = prob.select_entries(*[prob.DS, prob.GT, prob.pKO, prob.knockout, prob.G, prob.G_add_norm, prob.G_rec, prob.G_rec_norm_by_add, prob.G_rec_norm_by_rec])
    prob = prob.add_row_index()
    prob = prob.annotate_rows(
        locus = ht[prob.row_idx].locus,
        tmp_ref = ht[prob.row_idx].ref,
        tmp_alt = ht[prob.row_idx].alt,
        rsid=prob.gene_id
    )

    # conver to hail locus
    prob = prob.annotate_rows(
            locus=hl.parse_locus(prob.locus),
            alleles=[prob.tmp_ref, prob.tmp_alt]
    )

    # clean up and key appropiately
    prob = prob.key_rows_by(prob.locus, prob.alleles)
    prob = prob.drop(*['gene_id','tmp_ref','tmp_alt','row_idx'])

    # remove invariant sites
    prob = prob.annotate_rows(stdev = hl.agg.stats(prob.DS).stdev)
    prob = prob.filter_rows(prob.stdev > 0)

    # export genes with knockouts
    genes = prob.drop(prob.stats_add)
    genes = genes.drop(genes.stats_rec)
    genes = genes.drop(genes.stdev)
    genes.filter_entries(genes.DS > 0).entries().flatten().export(out_prefix + ".tsv.gz")

    # write out variants involved and vcf
    io.export_table(prob, out_prefix, out_type)

    # write matrix-table which contains dosages. VCF with only
    # dosages can't be re-read in HAIL, so we write a MatrixTable
    if out_prefix not in "mt":
        io.export_table(prob, out_prefix, "mt")



if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='Chromosome to be used')
    parser.add_argument('--in_prefix', default=None, help='Path to input')
    parser.add_argument('--in_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    # filtering options
    parser.add_argument('--maf_min', default=None, help='Select all variants with a maf greater than the indicated values')
    parser.add_argument('--maf_max', default=None, help='Select all variants with a maf less than the indicated value')
    parser.add_argument('--exclude', default=None, help='exclude variants by rsid and/or variant id')
    parser.add_argument('--use_loftee', default=False, action='store_true', help='use LOFTEE to distinghiush between high confidence PTVs')
    parser.add_argument('--exclude_singletons', default=False, action='store_true', help='Exclude singletons for analysis')
    parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='What categories should be subsetted to?')
    args = parser.parse_args()

    main(args)


