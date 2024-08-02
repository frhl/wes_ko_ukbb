#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import tables
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
    exclude = args.exclude
    csqs_category = args.csqs_category
   
    hail_init.hail_bmrc_init('logs/hail/append_vcf_rare.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
   
    mt = io.import_table(var_path, var_type, calc_info=False)
    
    if exclude:
        ht = hl.import_table(exclude, impute=True).key_by('varid')
        mt = mt.filter_rows(~hl.literal(set(ht.varid.collect())).contains(mt.varid))

    # Build variant annotation. Note: that we are building by variants,`
    # instead of gene, so that we get one variant per line. If we did this by
    # gene, in the scenario in which there are one variant affecting two different,
    # genes, then we would have the same variant twice (and these would be in perfect LD).
    
    # subset to current csqs category
    mt = mt.filter_rows(hl.literal(set(csqs_category)).contains(mt.consequence_category)) 

    # filter invariant sites
    mt = mt.annotate_entries(DS = hl.float(mt.GT.n_alt_alleles()))
  
    # export list of variants for later conditional analysis
    varid = hl.delimit([
        hl.str(mt.locus.contig),
        hl.str(mt.locus.position),
        mt.alleles[0],
        mt.alleles[1]],':')
    mt = mt.annotate_rows(rsid=varid)

    # get MAC
    mt = mt.annotate_rows(MAC=hl.min(hl.agg.call_stats(mt.GT, mt.alleles).AC))

    # create list of markers in data
    ht = mt.rows()
    ht = ht.select('rsid', 'consequence_category', 'MAC')
    ht.flatten().export(out_prefix + "_markers.txt.gz")

    # annotate dosage
    mt = mt.drop("GT")

    # load knockouts
    ko_mt = io.import_table(ko_path, ko_type, calc_info=False)
    ko_mt = tables.order_cols(ko_mt, mt)
    final = io.rbind_matrix_tables(mt, ko_mt )

    # Filter out invariant sites
    final = final.annotate_rows(stdev = hl.agg.stats(final.DS).stdev)
    final = final.filter_rows(final.stdev > 0)

    # return MatrixTable
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
    parser.add_argument('--exclude', default=None, help='exclude variants by rsid and/or variant id')
    parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='What categories should be subsetted to?')
    
    args = parser.parse_args()

    main(args)

