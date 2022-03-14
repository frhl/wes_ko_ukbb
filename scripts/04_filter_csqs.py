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
    sex = args.sex
    maf_max = (args.maf_max)
    maf_min = (args.maf_min)
    exclude = args.exclude
    use_loftee = args.use_loftee
    csqs_category = (args.csqs_category)

    # run parser
    hail_init.hail_bmrc_init('logs/hail/filter_csqs.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') 
    mt = io.import_table(input_path, input_type)
  
    if sex not in 'both':
        mt = samples.filter_to_sex(mt, sex)

    if maf_max and maf_min:
        mt = variants.filter_maf(mt, max_maf=float(maf_max),min_maf=float(maf_min))

    if exclude:
        ht = hl.import_table(exclude, impute = True).key_by('varid')
        mt = mt.filter_rows(~hl.literal(set(ht.varid.collect())).contains(mt.varid))

    # Build variant annotation
    mt = mt.explode_rows(mt.consequence.vep.worst_csq_by_gene_canonical)
    mt = mt.annotate_rows(
        consequence_category=ko.csqs_case_builder(
                worst_csq_expr=mt.consequence.vep.worst_csq_by_gene_canonical,
                use_loftee=use_loftee))    

    # subset to current csqs category
    category = "_".join(csqs_category)
    items = csqs_category
    print(items)
    mt = mt.filter_rows(hl.literal(set(items)).contains(mt.consequence_category)) 
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    
    parser.add_argument('--sex', default='both', help='Filter to sex (males or females)')
    parser.add_argument('--maf_min', default=None, help='Select all variants with a maf greater than the indicated values')
    parser.add_argument('--maf_max', default=None, help='Select all variants with a maf less than the indicated value')
    parser.add_argument('--exclude', default=None, help='exclude variants by rsid and/or variant id')
    parser.add_argument('--use_loftee', default=False, action='store_true', help='use LOFTEE to distinghiush between high confidence PTVs')
    parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='What categories should be subsetted to?')
    
    
    args = parser.parse_args()

    main(args)



