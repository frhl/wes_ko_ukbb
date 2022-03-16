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


def pseudo_variant(mt,  chrom):
    """ annotate psuedo variants which are required
    to export the as a Variant Call File
    """
    mt = mt.annotate_rows(
            locus=hl.parse_locus(str(chrom) + ':1'),
            alleles=hl.literal(['0', '1']))
    mt = mt.key_rows_by(mt.locus, mt.alleles)
    mt = mt.drop('gene_id')
    return(mt)



def main(args):
    
    chrom = args.chrom
    input_path = args.input_path
    input_type = args.input_type
    replicates = args.replicates
    seed = args.seed
    gene = args.gene
    checkpoint = args.checkpoint
    out_prefix = args.out_prefix
    out_type = args.out_type
    
    hail_init.hail_bmrc_init('logs/hail/array_permute.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
    mt = io.import_table(input_path, input_type, calc_info = False)
    
    if gene:
        mt = mt.filter_rows(mt.gene_id == gene)

    if checkpoint:
        mt = mt.checkpoint(checkpoint, overwrite = True)

    n = int(replicates)
    mts = list()
    for i in range(n):
        use_seed = int(seed) * i
        _mt = mt.annotate_entries(KO = hl.rand_bool(mt.pTKO, use_seed))
        _mt = _mt.annotate_entries(DS = ko.calc_prop_het_ko(_mt.KO, _mt.phased_het, _mt.unphased_het) * 2 )
        _mt = _mt.select_entries(_mt.DS)
        _mt = _mt.annotate_rows(k=i, seed=use_seed)
        _mt = _mt.annotate_rows(rsid = hl.delimit([_mt.gene_id, hl.str('-p'), hl.str(i)],'')) 
        mts.append(_mt)
    
    mt = hl.MatrixTable.union_rows(*mts)
    mt = pseudo_variant(mt, chrom)
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='')
    parser.add_argument('--input_path', default=None, help='')
    parser.add_argument('--input_type', default=None, help='')
    parser.add_argument('--gene', default=None, help='')
    parser.add_argument('--seed', default=None, help='Seed used for randomizing')
    parser.add_argument('--replicates', default=2, help='')
    parser.add_argument('--checkpoint', default=None, help='')
    parser.add_argument('--out_prefix', default=None, help='')
    parser.add_argument('--out_type', default=None, help='')
    args = parser.parse_args()

    main(args)



