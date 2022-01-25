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
    
    # parser
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    out_type   = args.out_type
    
    # variant filters
    chrom = int(args.chrom)
    af_max = (args.af_max)
    af_min = (args.af_min)
    maf_max = (args.maf_max)
    maf_min = (args.maf_min)
    missing = (args.missing)
    use_loftee = args.use_loftee
    csqs_category = (args.csqs_category)

    # sample filtering options
    sex = args.sex

    # run parser
    hail_init.hail_bmrc_init('logs/hail/knockout.log', 'GRCh38')
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
    csq_expr = mt.consequence.vep.worst_csq_by_gene_canonical 
    mt = mt.explode_rows(csq_expr)
    mt = mt.annotate_rows(
        consequence_cateogory=ko.csqs_case_builder(
                worst_csq_expr=csq_expr,
                use_loftee=use_loftee)
        )    

    # subset to current csqs category
    category = "_".join(csqs_category)
    items = csqs_category
    mt = mt.filter_rows(hl.literal(set(items)).contains(mt.consequence_category)) 
    print(f"chr{chrom}: evaluating '{category}' category")
    
    # convert to gene x sample matrix
    mt_gene=(mt.group_rows_by(csq_expr.gene_id)
            .aggregate(
                gts=hl.agg.collect(mt.GT),
                varid=hl.agg.collect(mt.varid),
                rsid=hl.agg.collect(mt.rsid)
                )
           ) 
     
    # determine genes that are knocked out
    mt_gene = ko.sum_gts_entries(mt_gene)
    expr_pko = ko.calc_prob_ko(mt_gene.hom_alt, mt_gene.phased, mt_gene.unphased)
    expr_ko = ko.annotate_knockout(mt_gene.hom_alt, expr_pko)
    mt_gene = mt_gene.annotate_entries(
            pKO = expr_pko,
            knockout = expr_ko
            )

    # convert to dosage and write vcf
    outfile_dosage = str(out_prefix) + "_" + str(category) + ".vcf.bgz"
    mt_vcf = mt_gene.annotate_entries(DS=mt_gene.pKO * 2)
    mt_vcf.select_entries(mt.DS).export(outfile_dosage)

    # Write to table 
    outfile_table = str(out_prefix) + "_" + str(category) + ".txt.gz"
    mt_gene.filter_entries(mt_gene.pKO >0).entries().flatten().export(outfile_table) 


    #outfile_ko_prob = str(out_prefix) + "_" + str(category) +'_ko_prob.tsv.bgz'
    #if export_ko_probability and not os.path.exists(outfile_ko_prob):
    #    mt_ko = analysis.gene_csqs_calc_pKO(mt1_subset, mt2_subset, 'dosage')
    #    mt_ko_entries = mt_ko.entries()
    #    mt_ko_entries = mt_ko_entries.filter(mt_ko_entries.pKO>0)
    #    mt_ko_entries.export(outfile_ko_prob)

    #outfile_ko_rsid = str(out_prefix) + "_" + str(category) + '_knockouts.tsv.bgz'
    #if export_ko_rsid and not os.path.exists(outfile_ko_rsid):
    #    mt_ko_rsid = analysis.gene_csqs_knockout_builder(mt1_subset)
    #    mt_ko_rsid.export(outfile_ko_rsid)

    #outfile_saige = str(out_prefix) + "_" + str(category) + "_ko"
    #if export_saige_vcf and not os.path.exists(outfile_saige):
    #    out = analysis.gene_csqs_calc_pKO_pseudoSNP(mt1_subset, mt2_subset, chrom)
    #    qc.export_table(out, out_prefix = out_prefix + "_" + category + "_ko", out_type = 'vcf')
    #    out.write(out_prefix + "_" + category + "_ko.mt")

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
    parser.add_argument('--sex', default=None, help='Filter to sex (males or females)')
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



