#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import io
from ko_utils import ko

def annotate_gene_knockouts(mt, gene_expr, chrom):
    """ annotate genes that are knockedout conditioned on current phase """
    
    genes = ko.aggr_phase_count_by_expr(mt, gene_expr)
    expr_pko = ko.calc_prob_ko(genes.hom_alt_n, genes.phased, genes.unphased)
    genes = genes.annotate_entries(pKO=expr_pko)

    prob = genes.annotate_entries(DS=genes.pKO * 2)
    prob = prob.select_entries(prob.DS)
    prob = prob.annotate_rows(
            locus=hl.parse_locus(str(chrom) + ':1'),
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
    mt = io.import_table(input_path, input_type)
    
    n = 2
    mts = list()
    for i in range(n):
        use_seed = int(seed) * i
        _mt = mt.transmute_entries(GT = ko.rand_flip_call(mt.GT, seed = use_seed))
        _gene_expr = _mt.consequence.vep.worst_csq_by_gene_canonical.gene_id 
        _mt = annotate_gene_knockouts(_mt, _gene_expr, chrom)
        _mt = _mt.annotate_rows(k=i, seed=use_seed)
        _mt = _mt.transmute_rows(rsid = hl.delimit([_mt.rsid, hl.str('-p'), hl.str(i)],'')) 
        mts.append(_mt)
        print(f"iter {i}")

    mt = hl.MatrixTable.union_rows(*mts)
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='')
    parser.add_argument('--input_path', default=None, help='')
    parser.add_argument('--input_type', default=None, help='')
    parser.add_argument('--permutations', default=None, help='')
    parser.add_argument('--max_permutations', default=None, help='')
    parser.add_argument('--seed', default=None, help='Seed used for randomizing')
    parser.add_argument('--n', default=2, help='')
    parser.add_argument('--out_prefix', default=None, help='')
    parser.add_argument('--out_type', default=None, help='')
    args = parser.parse_args()

    main(args)



