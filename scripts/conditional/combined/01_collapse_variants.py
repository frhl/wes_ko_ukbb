#!/usr/bin/env python3

import hail as hl
import argparse
import pandas as pd
import random
import string


from ukb_utils import hail_init
from ukb_utils import variants
from ko_utils import io
from ko_utils import ko


class SplitArgs(argparse.Action):
    r"""Method for splitting input csv into a list"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))


def get_tid(length=5):
    r"""method for getting random ID string for alleles"""
    return ''.join(random.choices(string.ascii_uppercase + string.ascii_lowercase, 
                                  k=length))

def main(args):

    mac_cutoff = args.mac_cutoff
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    out_type = args.out_type
    chrom = int(args.chrom)
    csqs_category = args.csqs_category

    # import phased/unphased data
    hail_init.hail_bmrc_init(log='logs/hail/knockout.log', default_reference='GRCh38', min_block_size=128)
    hl._set_flags(no_whole_stage_codegen='1')
    mt = io.import_table(input_path, input_type, calc_info = False)
    mt = mt.filter_rows(hl.literal(set(csqs_category)).contains(mt.consequence_category))
    
    # only collapse ultra rare variants
    if mac_cutoff:
        mt = mt.annotate_rows(MAC = variants.get_mac_expr(mt))
        mt = mt.filter_rows(mt.MAC <= int(mac_cutoff))
    
    # collapse to genes
    gene_expr = mt.consequence.vep.worst_csq_by_gene_canonical.gene_id
    genes = mt.group_rows_by(gene_expr).aggregate(DS=hl.agg.max(mt.GT.n_alt_alleles()))

    # setup sites and alleles
    rows = genes.count()[0]
    rng = range((rows)+1, (rows*2)+1)
    locus = [ "chr%s:%s" % (chrom, str(i+1)) for i in rng]
    ref = [ get_tid(4) for i in rng]
    alt = [ get_tid(4) for i in rng]

    # convert to HailTable
    df = pd.DataFrame({'locus':locus, 'ref':ref, 'alt':alt})
    ht = hl.Table.from_pandas(df)
    ht = ht.add_index()
    ht = ht.key_by('idx')

    # annotate knockout matrix
    prob = genes
    prob = prob.select_entries(prob.DS)
    prob = prob.add_row_index()
    prob = prob.annotate_rows(
        locus = ht[prob.row_idx].locus,
        tmp_ref = ht[prob.row_idx].ref, 
        tmp_alt = ht[prob.row_idx].alt,
        rsid=prob.gene_id.lower()
    )

    # conver to hail locus
    prob = prob.annotate_rows(
            locus=hl.parse_locus(prob.locus),
            alleles=[prob.tmp_ref, prob.tmp_alt]
    )

    # clean up and key appropiately
    prob = prob.key_rows_by(prob.locus, prob.alleles)
    prob = prob.drop(*['gene_id','tmp_ref','tmp_alt','row_idx'])

    # remove invariant sites
    prob = prob.annotate_rows(stdev = hl.agg.stats(prob.DS).stdev)
    prob = prob.filter_rows(prob.stdev > 0)

    # write matrix-table which contains dosages. VCF with only
    # dosages can't be re-read in HAIL, so we write a MatrixTable
    if out_type not in "mt":
        prob = prob.checkpoint(out_prefix + ".mt", overwrite=True)
    io.export_table(prob, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='Chromosome to be used') 
    parser.add_argument('--mac_cutoff', default=None, help='What mac_count should be used for collapsing') 
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    parser.add_argument('--csqs_category', default=None, action=SplitArgs, help='What categories should be subsetted to?')

    args = parser.parse_args()

    main(args)



