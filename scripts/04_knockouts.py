#!/usr/bin/env python3

import hail as hl
import argparse
import pandas
import os.path
import os

from gnomad.utils.vep import process_consequences
from ukb_utils import hail_init
from ukb_utils import genotypes
from ko_utils import qc
from ko_utils import analysis


class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))

def main(args):
    
    # parser
    input_phased_path = args.input_phased_path
    input_phased_type = args.input_phased_type
    input_unphased_path = args.input_unphased_path
    input_unphased_type = args.input_unphased_type
    out_prefix = args.out_prefix
    out_type   = args.out_type
    
    # variant filters
    chrom      = int(args.chrom)
    af_max     = (args.af_max)
    af_min     = (args.af_min)
    maf_max    = (args.maf_max)
    maf_min    = (args.maf_min)
    missing    = (args.missing)
    use_loftee = args.use_loftee
    csqs_category = (args.csqs_category)

    # sample filtering options
    sex = bool(args.sex)
    get_related = bool(args.get_related)
    get_unrelated = bool(args.get_unrelated)
    get_europeans = bool(args.get_europeans)

    # output options
    export_ko_probability = args.export_ko_probability
    export_saige_vcf = args.export_saige_vcf
    export_ko_rsid = args.export_ko_rsid

    # run parser
    hail_init.hail_bmrc_init('logs/hail/knockout.log', 'GRCh38')
    mt1 = qc.get_table(input_path=input_phased_path, input_type=input_phased_type) 
    mt2 = qc.get_table(input_path=input_unphased_path, input_type=input_unphased_type)
  
    if sex:
        mt1 = samples.filter_to_sex(mt, sex)
        mt2 = samples.filter_to_sex(mt, sex)

    if missing:
        mt1 = qc.filter_min_missing(mt1, float(missing))
        mt2 = qc.filter_min_missing(mt2, float(missing))
    
    if af_max:
        mt1 = qc.filter_max_af(mt1, float(af_max))
        mt2 = qc.filter_max_af(mt2, float(af_max))

    if af_min:
        mt1 = qf.fiter_min_af(mt1, float(af_min))
        mt2 = qf.fiter_min_af(mt2, float(af_min))

    if maf_max and maf_min:
        mt1 = qc.filter_maf(mt1, 
                max_maf = float(maf_max),
                min_maf = float(maf_min))
        mt2 = qc.filter_maf(mt2, 
                max_maf = float(maf_max),
                min_maf= float(maf_min))
        n1 = mt1.count()
        n2 = mt2.count()
        print(f"Filtering on maf_max={maf_max} and maf_min={maf_min}, resulting in {n1} and {n2}")
    
    # required before doing worst_csq_by_gene_canonical
    mt1 = mt1.explode_rows(mt1.consequence.vep.worst_csq_by_gene_canonical)
    mt2 = mt2.explode_rows(mt2.consequence.vep.worst_csq_by_gene_canonical)
    
    # get VEP annotation and add to rows
    by_gene_annotation1 = analysis.annotation_case_builder(mt1.consequence.vep.worst_csq_by_gene_canonical, use_loftee = use_loftee)
    by_gene_annotation2 = analysis.annotation_case_builder(mt2.consequence.vep.worst_csq_by_gene_canonical, use_loftee = use_loftee)
    mt1 = mt1.annotate_rows(consequence_category = by_gene_annotation1)    
    mt2 = mt2.annotate_rows(consequence_category = by_gene_annotation2)    

    # get current category
    category = "_".join(csqs_category)
    items = csqs_category

    print(f"chr{chrom}: evaluating '{category}' category")
    mt1_subset = mt1.filter_rows(hl.literal(set(items)).contains(mt1.consequence_category)) 
    mt2_subset = mt2.filter_rows(hl.literal(set(items)).contains(mt2.consequence_category))

    outfile_ko_prob = str(out_prefix) + "_" + str(category) +'_ko_prob.tsv.bgz'
    if export_ko_probability and not os.path.exists(outfile_ko_prob):
        mt_ko = analysis.gene_csqs_calc_pKO(mt1_subset, mt2_subset, 'dosage')
        mt_ko_entries = mt_ko.entries()
        mt_ko_entries = mt_ko_entries.filter(mt_ko_entries.pKO>0)
        mt_ko_entries.export(outfile_ko_prob)

    outfile_ko_rsid = str(out_prefix) + "_" + str(category) + '_knockouts.tsv.bgz'
    if export_ko_rsid and not os.path.exists(outfile_ko_rsid):
        mt_ko_rsid = analysis.gene_csqs_knockout_builder(mt1_subset)
        mt_ko_rsid.export(outfile_ko_rsid)

    outfile_saige = str(out_prefix) + "_" + str(category) + "_ko.vcf.bgz"
    if export_saige_vcf and not os.path.exists(outfile_saige):
        out = analysis.gene_csqs_calc_pKO_pseudoSNP(mt1_subset, mt2_subset, chrom)
        qc.export_table(out, out_prefix = out_prefix + "_" + category + "_ko", out_type = 'vcf')
        #undefined = out.aggregate_entries(hl.agg.sum(~hl.is_defined(out.DS)))
        #n = out.count()
        #print(f"chr{chrom}: undefined = {undefined}; variant/sample-count = {n}")


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    # required params
    parser.add_argument('--input_phased_path', default=None, help='Path to input')
    parser.add_argument('--input_phased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_unphased_path', default=None, help='Path to input that contains singletons')
    parser.add_argument('--input_unphased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--chrom', default=None, help='Chromosome to be used') 
    
    # variant and entry filtering
    parser.add_argument('--sex', default=None, help='Filter to sex')
    parser.add_argument('--maf_min', default=None, help='Select all variants with a maf greater than the indicated values')
    parser.add_argument('--maf_max', default=None, help='Select all variants with a maf less than the indicated value')
    parser.add_argument('--af_min', default=None, help='Select all variants with a AF greater than the indicated value')
    parser.add_argument('--af_max', default=None, help='Select all variants with a AF less than the indicated value')
    parser.add_argument('--missing', default=0.05, help='Filter variants by missingness threshold')
    parser.add_argument('--use_loftee', default=False, action='store_true', help='use LOFTEE to distinghiush between high confidence PTVs')
    parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='What categories should be subsetted to?')
 
    # sample filtering
    parser.add_argument('--get_related', action='store_true', help='Select all samples that are related')
    parser.add_argument('--get_unrelated', action='store_true', help='Select all samples that are unrelated') 
    parser.add_argument('--get_europeans', action='store_true', help='Filter to genetically confimed europeans.')
    
    # out 
    parser.add_argument('--export_ko_rsid', action='store_true', help='Exports the table with rsIDs involved in KOs.')
    parser.add_argument('--export_ko_probability', action='store_true', help='Exports the KO probability.')
    parser.add_argument('--export_saige_vcf', action='store_true', help='Export a "fake" VCF file that contains KO probabilities as DS field..')
    
    args = parser.parse_args()

    main(args)



