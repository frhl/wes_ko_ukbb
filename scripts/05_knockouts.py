#!/usr/bin/env python3

import hail as hl
import argparse
import pandas as pd

from ukb_utils import hail_init
from ukb_utils import samples
from ko_utils import variants
from ko_utils import io
from ko_utils import ko


class SplitArgs(argparse.Action):
    r"""Method for splitting input csv into a list"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))


def get_tid(length=5):
    r"""method for getting random ID string for alleles"""
    return ''.join(random.choices(string.ascii_uppercase + string.ascii_lowercase, 
                                  k=length))

def main(args):

    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    out_type = args.out_type
    chrom = int(args.chrom)
    only_vcf = args.only_vcf
    checkpoint = args.checkpoint
    aggr_method = args.aggr_method

    sex = args.sex
    maf_max = args.maf_max
    maf_min = args.maf_min
    exclude = args.exclude
    use_loftee = args.use_loftee
    export_all_gts = args.export_all_gts
    csqs_category = args.csqs_category

    # import phased/unphased data
    hail_init.hail_bmrc_init('logs/hail/knockout.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
    mt = io.import_table(input_path, input_type)

    # filtering
    if sex not in 'both':
        mt = samples.filter_to_sex(mt, sex)
    if maf_max and maf_min:
        mt = variants.filter_maf(mt, max_maf=float(maf_max), min_maf=float(maf_min))
    if exclude:
        ht = hl.import_table(exclude, impute=True).key_by('varid')
        mt = mt.filter_rows(~hl.literal(set(ht.varid.collect())).contains(mt.varid))

    # Build variant annotation
    mt = mt.explode_rows(mt.consequence.vep.worst_csq_by_gene_canonical)
    mt = mt.annotate_rows(
        consequence_category=ko.csqs_case_builder(
                worst_csq_expr=mt.consequence.vep.worst_csq_by_gene_canonical,
                use_loftee=use_loftee))

    # subset to current csqs category
    mt = mt.filter_rows(hl.literal(set(csqs_category)).contains(mt.consequence_category))

    # perform an aggregation based on "collect", which requires a lot
    # of memory but allows the variant ID to be returned as well alongside
    # with information of cis/trans-CHs and heterozygotes.
    gene_expr = mt.consequence.vep.worst_csq_by_gene_canonical.gene_id
    if aggr_method in "fast":
        genes = ko.aggr_phase_count_by_expr(mt, gene_expr)
        expr_pko = ko.calc_prob_ko(genes.hom_alt_n, genes.phased, genes.unphased)
        expr_ko = ko.annotate_knockout(genes.hom_alt_n, expr_pko)
    elif aggr_method in "collect":
        genes = ko.collect_phase_count_by_expr(mt, gene_expr)
        genes = ko.sum_gts_entries(genes)
        expr_pko = ko.calc_prob_ko(genes.hom_alt_n, genes.phased, genes.unphased)
        expr_ko = ko.annotate_knockout(genes.hom_alt_n, expr_pko, genes.phased)
    else:
        raise TypeError(str(aggr_method) + " is not allowed!")

    # calculate probability of being knocked out based on phased counts
    genes = genes.annotate_entries(pKO=expr_pko, knockout=expr_ko)

    # checkpoint for more effecient data use
    if checkpoint or aggr_method in "collect":
        genes = genes.checkpoint(out_prefix + "_checkpoint.mt", overwrite=True)

    # setup sites and alleles
    rows = genes.count()[0]
    locus = [ "chr%s:%s" % (chrom, str(i+1)) for i in range(rows)]
    ref = [ get_tid(4) for i in range(rows)]
    alt = [ get_tid(4) for i in range(rows)]

    # convert to HailTable
    df = pd.DataFrame({'locus':locus, 'ref':ref, 'alt':alt})
    ht = hl.Table.from_pandas(df)
    ht = ht.add_index()
    ht = ht.key_by('idx')

    # annotate knockout matrix
    prob = genes.annotate_entries(DS=genes.pKO * 2)
    prob = prob.select_entries(prob.DS)
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

    # write out variants involved and vcf
    io.export_table(prob, out_prefix, out_type)
    if not only_vcf:
        if export_all_gts:
            genes.filter_entries(hl.is_defined(genes.knockout)).entries().flatten().export(out_prefix + "_all.tsv.gz")
        else:
            genes.filter_entries(genes.pKO > 0).entries().flatten().export(out_prefix + ".tsv.gz")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--chrom', default=None, help='Chromosome to be used') 
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--only_vcf', default=False, action='store_true', help='Only return VCF (less memory required when running)')
    parser.add_argument('--checkpoint', default=False, action='store_true', help='Checkpoint gene-aggregation matrix to avoid Spark Memory overflow errors') 
    parser.add_argument('--aggr_method', default="collect", help='How should the CH matrix be generated?')
    # filtering options
    parser.add_argument('--sex', default='both', help='Filter to sex (males or females)')
    parser.add_argument('--maf_min', default=None, help='Select all variants with a maf greater than the indicated values')
    parser.add_argument('--maf_max', default=None, help='Select all variants with a maf less than the indicated value')
    parser.add_argument('--exclude', default=None, help='exclude variants by rsid and/or variant id')
    parser.add_argument('--use_loftee', default=False, action='store_true', help='use LOFTEE to distinghiush between high confidence PTVs')
    parser.add_argument('--export_all_gts', default=False, action='store_true', help='Exports a table of all csqs')
    parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='What categories should be subsetted to?')

    args = parser.parse_args()

    main(args)



