#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import tables
from ukb_utils import genotypes
from ukb_utils import samples
from ukb_utils import variants
from ko_utils import io

AUTOSOMES = list(map(str, range(1, 23)))

def main(args):
    
    ko_path = args.ko_path
    ko_type = args.ko_type
    markers = args.markers
    out_prefix = args.out_prefix
    out_type = args.out_type
    final_sample_list = args.final_sample_list
    
    hail_init.hail_bmrc_init('logs/hail/extract_markers.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
    print(markers)
    markers = markers.split(",") # format chr12:4214:A:T
    if len(markers) == 0:
        raise TypeError("no markers are present after splitting on ','")

    mt = genotypes.get_ukb_imputed_v3_bgen(AUTOSOMES)
    mt = variants.liftover(mt)
    mt = mt.annotate_rows(
            marker = hl.delimit([
                hl.str(mt.locus.contig), 
                hl.str(mt.locus.position),
                hl.str(mt.alleles[0]),
                hl.str(mt.alleles[1])
                ],':'))
    print(markers)
    print(mt.count())
    mt = mt.filter_rows(hl.literal(set(markers)).contains(mt.marker))
    mt = mt.annotate_entries(DS=hl.float64(mt.GT.n_alt_alleles()))
    mt = mt.select_entries(mt.DS)
    print(mt.count())
    
    if final_sample_list:
        ht_final_samples = hl.import_table(
                final_sample_list,
                no_header=True, key='f0',
                delimiter=',')
        mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key])) 

    # select only relevant rows
    mt = mt.select_rows(*[mt.rsid, mt.varid])
    
    # extract indiviudal GT entries
    mt.entries().flatten().export(out_prefix + ".tsv.gz")

    # export variants to be used for conditional analysis
    io.export_table(mt, out_prefix, out_type)
    





if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ko_path', default=None, help='')
    parser.add_argument('--ko_type', default=None, help='')
    parser.add_argument('--markers', default=None, help='')
    parser.add_argument('--out_prefix', default=None, help='')
    parser.add_argument('--out_type', default=None, help='')
    parser.add_argument('--final_sample_list', default=None, help='')
    args = parser.parse_args()

    main(args)



