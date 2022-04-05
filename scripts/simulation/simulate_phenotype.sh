#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ko_utils import ko
from ukb_utils import hail_init

class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))

def main(args):

    in_prefix = args.in_prefix
    in_type = args.in_type
    out_prefix = args.out_prefix
    out_type = args.out_type
    chrom = args.chrom
    csqs_category = args.csqs_category
    h2 = args.h2
    pi = args.pi
    simulations = args.simulations
    export_single_markers = args.export_single_markers
    seed = args.seed
    
    # set parameters
    pi = float(pi) if pi is not None else None
    h2 = float(h2) if h2 is not None else None
 
    # import table
    hail_init.hail_bmrc_init('logs/hail/absence_of_effect.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
    hl.set_global_seed(int(seed))
    mt = io.import_table(in_prefix, in_type)

    # annotate with variant consequence
    mt = mt.explode_rows(mt.consequence.vep.worst_csq_by_gene_canonical)
    mt = mt.annotate_rows(
    consequence_category=ko.csqs_case_builder(
            worst_csq_expr=mt.consequence.vep.worst_csq_by_gene_canonical,
            use_loftee=True))
    # simulate betas
    mt = hl.experimental.ldscsim.make_betas(mt, float(h2), pi = pi)[0]
    mt = mt.filter_rows(hl.agg.stats(mt.GT.n_alt_alleles()).stdev>0)
    mt = hl.experimental.ldscsim.normalize_genotypes(mt.GT)  

    # add noise/environemental component
    mt = mt.annotate_cols(y_no_noise=hl.agg.sum(mt.beta * mt['norm_gt']))

    # simulate x times
    for i in range(int(simulations)):
        col = "y" + str(i)
        y = mt.y_no_noise + hl.rand_norm(0, hl.sqrt(1-h2))
        mt = mt.annotate_cols(y = y)
        mt = mt.rename({'y': col})

    # export simulated phenotypes
    ht = mt.cols()
    ht.flatten().export(out_prefix + ".tsv.gz")
    
    if export_single_markers:
        io.export_table(mt, out_prefix + "_markers", out_type)

    # collapse to gene knockout
    mt = mt.filter_rows(hl.literal(set(csqs_category)).contains(mt.consequence_category))
    gene_expr = mt.consequence.vep.worst_csq_by_gene_canonical.gene_id
    genes = ko.aggr_phase_count_by_expr(mt, gene_expr)
    
    # calculate p(KO)
    expr_pko = ko.calc_prob_ko(genes.hom_alt_n, genes.phased, genes.unphased)
    expr_ko = ko.annotate_knockout(genes.hom_alt_n, expr_pko)
    genes = genes.annotate_entries(
            pKO=expr_pko,
            knockout=expr_ko)

    # annotate dosage (used by SAIGE)
    prob = genes.annotate_entries(DS=genes.pKO * 2)
    prob = prob.select_entries(prob.DS)
    prob = prob.annotate_rows(
            locus=hl.parse_locus('chr' + str(chrom) + ':1'),
            alleles=hl.literal(['0', '1']),
            rsid=prob.gene_id)
    prob = prob.key_rows_by(prob.locus, prob.alleles)
    prob = prob.drop('gene_id')
    io.export_table(prob, out_prefix, out_type)

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='chromosome')
    parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='comma sepearted field')
    parser.add_argument('--h2', default=None, help='Heritability for phenotype simulated')
    parser.add_argument('--pi', default=None, help='Probability of variant being causal')
    parser.add_argument('--simulations', default=1, help='simulations to be dobe')
    parser.add_argument('--seed', default=None, help='seed for random simulations')
    parser.add_argument('--in_prefix', default=None, help='Path prefix for input dataset')
    parser.add_argument('--in_type', default=None, help='Either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Either "mt", "vcf" or "plink"')
    parser.add_argument('--export_single_markers', action='store_true', help='Either "mt", "vcf" or "plink"')
    args = parser.parse_args()

    main(args)


