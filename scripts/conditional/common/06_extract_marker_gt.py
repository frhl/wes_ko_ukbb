#!/usr/bin/env python3

import hail as hl
import argparse

import pandas as pd
from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import variants
from ko_utils import io

AUTOSOMES = list(map(str, range(1, 23)))

def main(args):
    
    markers = args.markers
    out_prefix = args.out_prefix
    out_type = args.out_type
    final_sample_list = args.final_sample_list
    
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
   
    # assume input is grch38 and translate to grch37 
    from_build = 'GRCh38'
    to_build = 'GRCh37'
    liftover_path = variants.get_liftover_chain_path(from_build,to_build)
    rg_from = hl.get_reference(from_build)  
    rg_to = hl.get_reference(to_build)
    if not rg_from.has_liftover(rg_to):
        rg_from.add_liftover(liftover_path, rg_to)   

    # liftover to grch37
    ht = ht.annotate(new_locus=hl.liftover(ht.locus, to_build, include_strand=True),old_locus=ht.locus)
    ht = ht.filter(hl.is_defined(ht.new_locus))
    ht = ht.annotate(new_alleles=hl.if_else(ht.new_locus.is_negative_strand, [variants.flip_base(ht.alleles[0]), variants.flip_base(ht.alleles[1])],ht.alleles))
    ht = ht.key_by(locus=ht.new_locus.result, alleles=ht.new_alleles)
   
    # get chromosomes aad subset imputed data to varaints in hailTable
    chroms = set(ht.locus.contig.collect())
    mt = genotypes.get_ukb_imputed_v3_bgen(chroms)
    mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))
    mt = mt.annotate_rows(grch38 = ht[mt.row_key].grch38)
    mt = mt.key_rows_by(locus=mt.grch38.locus, alleles=mt.grch38.alleles)     

    # We are interested in dosages
    mt = mt.annotate_entries(DS=hl.float64(mt.GT.n_alt_alleles()))
    mt = mt.select_entries(*[mt.DS, mt.GT])
    mt = mt.drop(mt.grch38)

    # fitler to releveant samples
    if final_sample_list:
        ht_final_samples = hl.import_table(
                final_sample_list,
                no_header=True, key='f0',
                delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))

    # export variants to be used for conditional analysis
    io.export_table(mt, out_prefix, "mt")
    mt = io.import_table(out_prefix + ".mt", "mt")

    # write VCF
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
    parser.add_argument('--final_sample_list', default=None, help='')
    args = parser.parse_args()

    main(args)



