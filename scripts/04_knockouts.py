#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import samples
from ko_utils import variants
from ko_utils import io
from ko_utils import ko


class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))

def main(args):
    
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    out_type = args.out_type
    
    chrom = int(args.chrom)
    maf_max = (args.maf_max)
    maf_min = (args.maf_min)
    exclude = args.exclude
    use_loftee = args.use_loftee
    only_vcf = args.only_vcf
    csqs_category = (args.csqs_category)
    randomize_phase = args.randomize_phase
    seed = args.seed
    aggr_method = args.aggr_method
    checkpoint = args.checkpoint
    sex = args.sex

    # run parser
    hail_init.hail_bmrc_init('logs/hail/knockout.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') 
    mt = io.import_table(input_path, input_type)
  
    if sex not in 'both':
        mt = samples.filter_to_sex(mt, sex)

    if maf_max and maf_min:
        mt = variants.filter_maf(mt, max_maf=float(maf_max),min_maf=float(maf_min))
   
    if exclude:
        exlcude = exclude.split(',').strip()
        rsid_expr = hl.literal(set(exclude)).contains(mt.rsid)
        varid_expr = hl.literal(set(exclude)).contains(mt.varid)
        mt = mt.filter_rows(~(rsid_expr | varid_expr))

    if randomize_phase:
        hetz_before = ko.aggr_count_calls(mt)
        mt = mt.transmute_entries(GT=ko.rand_flip_call(mt.GT, seed=int(seed)))
        hetz_after = ko.aggr_count_calls(mt)
        print(f"Phased hetz before {hetz_before} and after {hetz_after}")

    # Build variant annotation
    mt = mt.explode_rows(mt.consequence.vep.worst_csq_by_gene_canonical)
    mt = mt.annotate_rows(
        consequence_category=ko.csqs_case_builder(
                worst_csq_expr=mt.consequence.vep.worst_csq_by_gene_canonical,
                use_loftee=use_loftee))    

    # subset to current csqs category
    category = "_".join(csqs_category)
    items = csqs_category
    mt = mt.filter_rows(hl.literal(set(items)).contains(mt.consequence_category)) 

    # perform an aggregation based on "collect", which requires a lot
    # of memory but allows the variant ID to be returned as well?
    gene_expr = mt.consequence.vep.worst_csq_by_gene_canonical.gene_id
    if aggr_method in "fast":
        print("[Note]: Using fast aggregation..")
        genes = ko.aggr_phase_count_by_expr(mt, gene_expr)
    elif aggr_method in "collect":
        genes = ko.collect_phase_count_by_expr(mt, gene_expr) 
        genes = ko.sum_gts_entries(genes)
    else:
        raise TypeError(str(aggr_method) + " is not allowed!") 

    # calculate probability of being knocked out based on phased counts
    expr_pko = ko.calc_prob_ko(genes.hom_alt_n, genes.phased, genes.unphased)
    expr_ko = ko.annotate_knockout(genes.hom_alt_n, expr_pko)
    genes = genes.annotate_entries(
            pKO=expr_pko,
            knockout=expr_ko)

    # checkpoint for more effecient data use
    if checkpoint or aggr_method in "collect":
        genes = genes.checkpoint(out_prefix + "_checkpoint.mt", overwrite=True)

    # get probabilistic matrix for knockouts
    prob = genes.annotate_entries(DS=genes.pKO * 2)
    prob = prob.select_entries(prob.DS)
    prob = prob.annotate_rows(
            locus=hl.parse_locus('chr' + str(chrom) + ':1'),
            alleles=hl.literal(['0', '1']),
            rsid=prob.gene_id)
    prob = prob.key_rows_by(prob.locus, prob.alleles)
    prob = prob.drop('gene_id')
    
    # write out 
    csq_prefix = str(out_prefix) + "_" + str(category)
    io.export_table(prob, csq_prefix, out_type)
    if not only_vcf:
        genes.filter_entries(genes.pKO > 0).entries().flatten().export(csq_prefix + ".tsv.gz") 

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--chrom', default=None, help='Chromosome to be used') 
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--aggr_method', default="collect", help='How should the CH matrix be generated?')
    
    parser.add_argument('--sex', default='both', help='Filter to sex (males or females)')
    parser.add_argument('--maf_min', default=None, help='Select all variants with a maf greater than the indicated values')
    parser.add_argument('--maf_max', default=None, help='Select all variants with a maf less than the indicated value')
    parser.add_argument('--exclude', default=None, help='exclude variants by rsid and/or variant id')
    parser.add_argument('--use_loftee', default=False, action='store_true', help='use LOFTEE to distinghiush between high confidence PTVs')
    parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='What categories should be subsetted to?')
    
    parser.add_argument('--only_vcf', default=False, action='store_true', help='Only return VCF (less memory required when running)')
    parser.add_argument('--checkpoint', default=False, action='store_true', help='Checkpoint gene-aggregation matrix to avoid Spark Memory overflow errors') 
    parser.add_argument('--randomize_phase', default=None, action='store_true', help='Randomize phased calls?')
    parser.add_argument('--seed', default=None, help='Seed used for randomizing')
    
    args = parser.parse_args()

    main(args)



