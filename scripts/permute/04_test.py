#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import io
from ko_utils import ko

def annotate_gene_knockouts(mt, gene_expr):
    """ annotate genes that are knockedout conditioned on current phase """
    
    genes = ko.aggr_phase_count_by_expr(mt, gene_expr)
    expr_pko = ko.calc_prob_ko(genes.hom_alt_n, genes.phased, genes.unphased)
    genes = genes.annotate_entries(pKO=expr_pko)

    prob = genes.annotate_entries(DS=genes.pKO * 2)
    prob = prob.select_entries(prob.DS)
    prob = prob.annotate_rows(
            locus=hl.parse_locus('chr' + str(chrom) + ':1'),
            alleles=hl.literal(['0', '1']),
            rsid=prob.gene_id)
    prob = prob.key_rows_by(prob.locus, prob.alleles)
    prob = prob.drop('gene_id')
    return(prob)


def main(args):
    
    chrom = args.chrom
    input_path = args.input_path
    input_type = args.input_type
    permutations = args.permutations
    max_permutations = args.max_permutations
    seed = args.seed
    out_prefix = args.out_prefix
    out_type = args.out_type
    
    hail_init.hail_bmrc_init('logs/hail/array_permute.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')

    # get table of variants involed in knockouts
    ht = hl.import_table(
        permutations,
        force = True,
        impute = True
    )

    # get variants in the speceific chromosome
    ht = ht.filter(ht.CHR == chrom)
    variants = ", ".join(ht.varid.collect()).split(', ') 

    # how many permutations required?
    n_cur = int(ht.aggregate(hl.agg.max(ht.permut)))
    n = min(n_cur, int(max_permutations))
    print(f"{n_cur} permutations required. Using {n}.")

    # get matrix table with variants
    mt = io.import_table(input_path, input_type)
    
    # subet to damaging category
    category = "_".join(csqs_category)
    items = csqs_category
    mt = mt.filter_rows(hl.literal(set(items)).contains(mt.consequence_category))
    mt = mt.checkpoint(out_prefix + "_checkpoint.mt", overwrite=True)
    gene_expr = mt.consequence.vep.worst_csq_by_gene_canonical.gene_id
    
       
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='')
    parser.add_argument('--input_path', default=None, help='')
    parser.add_argument('--input_type', default=None, help='')
    parser.add_argument('--permutations', default=None, help='')
    parser.add_argument('--max_permutations', default=None, help='')
    parser.add_argument('--seed', default=None, help='Seed used for randomizing')
    parser.add_argument('--out_prefix', default=None, help='')
    parser.add_argument('--out_type', default=None, help='')
    args = parser.parse_args()

    main(args)



