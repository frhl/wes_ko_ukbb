#!/usr/bin/env python3

import hail as hl
import argparse
import pandas
import os

from gnomad.utils.vep import process_consequences
from ukb_utils import hail_init
from ukb_utils import genotypes
from ko_utils import qc
from ko_utils import analysis

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
    hwe        = (args.hwe)
    missing    = (args.missing)
    min_dp     = (args.min_dp)
    min_gq     = (args.min_gq)
    
    # sample filtering options
    get_related = bool(args.get_related)
    get_unrelated = bool(args.get_unrelated)
    get_europeans = bool(args.get_europeans)

    # output options
    vep_filter = args.vep_filter
    export_ko_dosage_matrix = args.export_ko_dosage_matrix
    export_burden = args.export_burden
    export_ko_probability = args.export_ko_probability
    export_saige_vcf = args.export_saige_vcf
    export_ko_rsid = args.export_ko_rsid

    # run parser
    hail_init.hail_bmrc_init('logs/hail/knockout.log', 'GRCh38')
    mt1 = qc.get_table(input_path=input_phased_path, input_type=input_phased_type) 
    mt2 = qc.get_table(input_path=input_unphased_path, input_type=input_unphased_type)
    
    # add tmp rsid
    mt1 = mt1.annotate_rows(rsid = mt1.snpid)
    mt2 = mt2.annotate_rows(rsid = mt2.snpid)

    if get_related and not get_unrelated:
        mt1 = qc.filter_to_unrelated(mt1, get_related = True)
        mt2 = qc.filter_to_unrelated(mt2, get_related = True)

    if get_unrelated and not get_related:
        mt1 = qc.filter_to_unrelated(mt1, get_related = False)
        mt2 = qc.filter_to_unrelated(mt2, get_related = False)
        
    if get_europeans: 
        mt1 = qc.filter_to_european(mt1, genetically_european = True)
        mt2 = qc.filter_to_european(mt2, genetically_european = True)

    if missing:
        mt1 = qc.filter_min_missing(mt1, float(missing))
        mt2 = qc.filter_min_missing(mt2, float(missing))
    
    if af_max:
        mt1 = qc.filter_max_af(mt1, float(af_max))

    if af_min:
        mt1 = qf.fiter_min_af(mt1, float(af_min))

    if maf_max:
        mt1 = qc.filter_max_maf(mt1, float(maf_max))

    if maf_min:
        mt1 = qc.filter_min_maf(mt1, float(maf_min))

    if hwe:
        mt1 = qc.filter_hwe(mt1, float(hwe))

    if min_dp:
        mt1 = mt1
        #mt1 = mt1.filter_entries(mt1.DP >= min_dp) 
        #mt2 = mt2.filter_entries(mt2.DP >= min_dp)

    if min_gq:
        mt1 = mt1
        #mt1 = mt1.filter_entries(mt1.GQ >= min_gq) 
        #mt1 = mt1.filter_entries(mt1.GQ >= min_gq) 

    # Annotate burden variant category
    mt1 = mt1.explode_rows(mt1.vep.worst_csq_by_gene_canonical)
    mt1 = analysis.variant_csqs_category_builder(mt1)
    mt2 = mt2.explode_rows(mt2.vep.worst_csq_by_gene_canonical)
    mt2 = analysis.variant_csqs_category_builder(mt2)

    if vep_filter: 
        mt1 = analysis.filter_vep(mt1, 'consequence_category', vep_filter)
        mt2 = analysis.filter_vep(mt2, 'consequence_category', vep_filter) 

    #### get stats

    if export_burden:

        # Count burden per gene per individual
        mt1_cat = analysis.gene_burden_category_annotations_per_sample(mt1)
        mt2_cat = analysis.gene_burden_category_annotations_per_sample(mt2)

        # combine singleton table and full table
        res = mt1_cat.annotate_entries(singletons = mt2_cat[(mt1_cat.gene_id, mt1_cat.consequence_category), mt1_cat.s].n)
        res = res.annotate_entries(singletons = hl.if_else(hl.is_missing(res.singletons),0,res.singletons))
        res = res.annotate_entries(total = res.n + res.singletons)
        res = res.entries()
        res = res.filter(res.total > 0)

        # export data
        res.export(out_prefix + '_burden.tsv.gz')

    if export_ko_probability:

        # determine probability of being a ko
        mt_ko = analysis.gene_csqs_calc_pKO(mt1, mt2, 'dosage')
        mt_ko_entries = mt_ko.entries()
        mt_ko_entries = mt_ko_entries.filter(mt_ko_entries.pKO>0)

        # export data
        mt_ko_entries.export(out_prefix + '_ko_prob.tsv.gz')

    if export_ko_rsid:
        mt_ko_rsid = analysis.gene_csqs_knockout_builder(mt1)
        mt_ko_rsid.export(out_prefix + '_knockouts.tsv.gz')

    if export_saige_vcf:
        out = analysis.gene_csqs_calc_pKO_pseudoSNP(mt1, mt2, chrom)
        qc.export_table(out, out_prefix = out_prefix + "_ko", out_type = 'vcf')

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
    parser.add_argument('--maf_min', default=None, help='Select all variants with a maf greater than the indicated values')
    parser.add_argument('--maf_max', default=None, help='Select all variants with a maf less than the indicated value')
    parser.add_argument('--af_min', default=None, help='Select all variants with a AF greater than the indicated value')
    parser.add_argument('--af_max', default=None, help='Select all variants with a AF less than the indicated value')
    parser.add_argument('--hwe', default=None, help='Filter variants by HWE threshold')
    parser.add_argument('--missing', default=0.05, help='Filter variants by missingness threshold')
    parser.add_argument('--min_dp', default=20, help='filter variants by minimum sequencing depth (DP)')
    parser.add_argument('--min_gq', default=48, help='filter variants by genotype quality (GQ)')
    
    # sample filtering
    parser.add_argument('--get_related', action='store_true', help='Select all samples that are related')
    parser.add_argument('--get_unrelated', action='store_true', help='Select all samples that are unrelated') 
    parser.add_argument('--get_europeans', action='store_true', help='Filter to genetically confimed europeans.')
    
    # out 
    parser.add_argument('--export_ko_rsid', action='store_true', help='Exports the table with rsIDs involved in KOs.')
    parser.add_argument('--export_ko_probability', action='store_true', help='Exports the KO probability.')
    parser.add_argument('--export_burden', action='store_true', help='Export burden variant count by gene and and individuals.')
    parser.add_argument('--export_saige_vcf', action='store_true', help='Export a "fake" VCF file that contains KO probabilities as DS field..')
    parser.add_argument('--vep_filter', nargs='+', help='Filter consequence_category by mutations e.g., "damaging_missense" or "ptv"')
    parser.add_argument('--export_ko_dosage_matrix', action='store_true', help='Generate a gene x sample matrix with KO status')
    
    args = parser.parse_args()

    main(args)



