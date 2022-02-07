#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import samples
from ko_utils import variants as va
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
    af_max = (args.af_max)
    af_min = (args.af_min)
    maf_max = (args.maf_max)
    maf_min = (args.maf_min)
    missing = (args.missing)
    use_loftee = args.use_loftee
    csqs_category = (args.csqs_category)

    sex = args.sex

    # run parser
    hail_init.hail_bmrc_init('logs/hail/knockout.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') 
    mt = io.import_table(
            input_path=input_path,
            input_type=input_type
            ) 
  
    if sex:
        if sex not in 'both':
            mt = samples.filter_to_sex(mt, sex)

    if missing:
        mt = va.filter_min_missing(mt, float(missing))
    
    if af_max:
        mt = va.filter_max_af(mt, float(af_max))

    if af_min:
        mt = va.fiter_min_af(mt, float(af_min))

    if maf_max and maf_min:
        mt = va.filter_maf(mt, 
                max_maf = float(maf_max),
                min_maf = float(maf_min))
        n = mt.count()
        print(f"${n} samples remaining after filtering on maf_max={maf_max} and maf_min={maf_min}")
  
    # Build variant annotation
    mt = mt.explode_rows(mt.consequence.vep.worst_csq_by_gene_canonical)
    mt = mt.annotate_rows(
        consequence_category=ko.csqs_case_builder(
                worst_csq_expr=mt.consequence.vep.worst_csq_by_gene_canonical,
                use_loftee=use_loftee)
        )    

    # subset to current csqs category
    
    category = "_".join(csqs_category)
    items = csqs_category
    mt = mt.filter_rows(hl.literal(set(items)).contains(mt.consequence_category)) 
    print(f"chr{chrom}: evaluating '{category}' category")
    
    # convert to gene x sample matrix
    mt_gene = (mt.group_rows_by(mt.consequence.vep.worst_csq_by_gene_canonical.gene_id)
            .aggregate(
                gts=hl.agg.collect(mt.GT),
                varid=hl.agg.collect(mt.varid)
                #rsid=hl.agg.collect(mt.rsid)
                )
           )

    # determine genes that are knocked out
    mt_gene = ko.sum_gts_entries(mt_gene)
    expr_pko = ko.calc_prob_ko(mt_gene.hom_alt, mt_gene.phased, mt_gene.unphased)
    expr_ko = ko.annotate_knockout(mt_gene.hom_alt, expr_pko)
    mt_gene = mt_gene.annotate_entries(
            pKO = expr_pko,
            knockout = expr_ko)

    # convert to dosage and write vcf
    csq_prefix = str(out_prefix) + "_" + str(category)
    
    mt_vcf = mt_gene.annotate_entries(DS=mt_gene.pKO * 2)
    mt_vcf = mt_vcf.select_entries(mt_vcf.DS)
    mt_vcf = mt_vcf.annotate_rows(
            locus=hl.parse_locus('chr' + str(chrom) + ':1'),
            alleles=hl.literal(['X', 'Y']),
            rsid=mt_vcf.gene_id)

    mt_vcf = mt_vcf.key_rows_by(mt_vcf.locus, mt_vcf.alleles)
    mt_vcf = mt_vcf.drop('gene_id')
    hl.export_vcf(mt_vcf, csq_prefix + '.vcf.bgz')
    io.export_table(mt_vcf, csq_prefix, out_type)

    # Write to table
    mt_gene.filter_entries(mt_gene.pKO > 0).entries().flatten().export(csq_prefix + "tsv.gz") 

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--chrom', default=None, help='Chromosome to be used') 
    
    parser.add_argument('--sex', default=None, help='Filter to sex (males or females)')
    parser.add_argument('--maf_min', default=None, help='Select all variants with a maf greater than the indicated values')
    parser.add_argument('--maf_max', default=None, help='Select all variants with a maf less than the indicated value')
    parser.add_argument('--af_min', default=None, help='Select all variants with a AF greater than the indicated value')
    parser.add_argument('--af_max', default=None, help='Select all variants with a AF less than the indicated value')
    parser.add_argument('--missing', default=0.05, help='Filter variants by missingness threshold')
    parser.add_argument('--use_loftee', default=False, action='store_true', help='use LOFTEE to distinghiush between high confidence PTVs')
    parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='What categories should be subsetted to?')
 
    parser.add_argument('--get_related', action='store_true', help='Select all samples that are related')
    parser.add_argument('--get_unrelated', action='store_true', help='Select all samples that are unrelated') 
    parser.add_argument('--get_europeans', action='store_true', help='Filter to genetically confimed europeans.')
    
    parser.add_argument('--export_ko_rsid', action='store_true', help='Exports the table with rsIDs involved in KOs.')
    parser.add_argument('--export_ko_probability', action='store_true', help='Exports the KO probability.')
    parser.add_argument('--export_saige_vcf', action='store_true', help='Export a "fake" VCF file that contains KO probabilities as DS field..')
    
    args = parser.parse_args()

    main(args)



