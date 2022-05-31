#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import tables
from ukb_utils import samples
from ko_utils import variants
from ko_utils import io
from ko_utils import ko

class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))

def main(args):
    
    ko_path = args.ko_path
    var_path = args.var_path
    ko_type = args.ko_type
    var_type = args.var_type
    out_prefix = args.out_prefix
    out_type = args.out_type
    sex = args.sex
    maf_max = args.maf_max
    maf_min = args.maf_min
    exclude = args.exclude
    use_loftee = args.use_loftee
    csqs_category = args.csqs_category
   
    hail_init.hail_bmrc_init('logs/hail/append_vcf_rare.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
   
    mt = io.import_table(var_path, var_type, calc_info=False)
    
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
    
    # export list of variants for later conditional analysis
    varid = hl.delimit([
        hl.str(mt.locus.contig),
        hl.str(mt.locus.position),
        mt.alleles[0],
        mt.alleles[1]],':')
    mt = mt.annotate_rows(rsid=varid)
    ht = mt.rows()
    ht = ht.select('rsid', 'consequence_category')
    ht.flatten().export(out_prefix + "_markers.txt.gz")

    # annotate dosage
    mt = mt.annotate_entries(DS = hl.float(mt.GT.n_alt_alleles()))
    mt = mt.drop("GT")

    # load knockouts
    ko_mt = io.import_table(ko_path, ko_type, calc_info=False)
    ko_mt = tables.order_cols(ko_mt, mt)
    final = io.rbind_matrix_tables(mt, ko_mt )
   
    if out_type not in "mt":
        io.export_table(final, out_prefix, "mt")
    io.export_table(final, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--var_path', default=None, help='')
    parser.add_argument('--var_type', default=None, help='')
    parser.add_argument('--ko_path', default=None, help='')
    parser.add_argument('--ko_type', default=None, help='')
    parser.add_argument('--out_prefix', default=None, help='')
    parser.add_argument('--out_type', default=None, help='')
    
    parser.add_argument('--sex', default='both', help='Filter to sex (males or females)')
    parser.add_argument('--maf_min', default=None, help='Select all variants with a maf greater than the indicated values')
    parser.add_argument('--maf_max', default=None, help='Select all variants with a maf less than the indicated value')
    parser.add_argument('--exclude', default=None, help='exclude variants by rsid and/or variant id')
    parser.add_argument('--use_loftee', default=False, action='store_true', help='use LOFTEE to distinghiush between high confidence PTVs')
    parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='What categories should be subsetted to?')
    
    args = parser.parse_args()

    main(args)

