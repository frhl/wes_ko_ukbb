#!/usr/bin/env python3

import hail as hl
import argparse

import pandas as pd
from ukb_utils import hail_init
from ko_utils import io


def main(args):
    
    markers = args.markers
    out_prefix = args.out_prefix
    out_type = args.out_type
    
    hail_init.hail_bmrc_init('logs/hail/extract_markers.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
    markers = markers.split(",") # format chr12:4214:A:T
    markers = list(filter(None, markers)) 
    print(markers)
    if len(markers) == 0:
        raise TypeError("no markers are present after splitting on ','")

    # read in table of markers
    d = {'input': markers}
    ht = hl.Table.from_pandas(pd.DataFrame(data=d))
    ht = ht.annotate(var = hl.parse_variant(ht.input))
    ht = ht.annotate(locus = ht.var.locus, alleles = ht.var.alleles)
    ht = ht.drop("var")

    # annotate marker / variant category
    ht = ht.annotate(
        csqs = "common",
        marker = hl.delimit([
            hl.str(ht.locus.contig),
            hl.str(ht.locus.position),
            hl.str(ht.alleles[0]),
            hl.str(ht.alleles[1])
            ],':'))

    # keep grch38 information
    ht = ht.annotate(
        grch38 = hl.struct(
            locus = ht.locus,
            alleles = ht.alleles,
            marker = ht.marker
        )
    )
    
    # get chromosomes to import
    ht = ht.key_by(locus=ht.locus, alleles=ht.alleles)
    chromosomes = set(ht.locus.contig.collect())
    
    mts = list()
    for chrom in chromosomes:
        path = f"data/unphased/imputed/common_append_missing/ukb_imp_200k_common_append_missing_{chrom}.mt"
        mt = hl.read_matrix_table(path)
        mts.append(mt)
    
    # combine and subset
    mt = mts[0].union_rows(*mts[1:]) 
    mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))
    
    # We are interested in dosages for SAIGE
    mt = mt.annotate_entries(DS=hl.float64(mt.GT.n_alt_alleles()))
    mt = mt.select_entries(*[mt.DS, mt.GT])

    # export variants to be used for conditional analysis
    io.export_table(mt, out_prefix, "mt")
    if out_prefix not in "mt":
        io.export_table(mt, out_prefix, out_type)
    
    # extract indiviudal GT entries
    mt = mt.select_rows(*[mt.rsid, mt.varid])
    mt.entries().flatten().export(out_prefix + ".tsv.gz")

    # get ordered table of variants included
    ht = ht.flatten()
    ht = ht.select(*[ht.locus, ht.alleles, ht.marker, ht.csqs])
    ht = ht.rename({"marker" : "rsid"})
    ht.flatten().export(out_prefix + "_rows.txt.gz")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--markers', default=None, help='')
    parser.add_argument('--out_prefix', default=None, help='')
    parser.add_argument('--out_type', default=None, help='')
    args = parser.parse_args()

    main(args)



