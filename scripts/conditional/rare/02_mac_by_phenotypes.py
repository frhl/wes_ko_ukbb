#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import io

class SplitArgsByComma(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))

def main(args):
    
    input_path = args.input_path
    input_type = args.input_type
    pheno_path = args.pheno_path
    phenotypes = args.phenotypes
    out_prefix = args.out_prefix
    out_type = args.out_type

    hail_init.hail_bmrc_init('logs/hail/append_vcf_rare.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')

    mt = io.import_table(input_path, input_type, calc_info = False)
    ht = hl.import_table(pheno_path,
             types={'eid': hl.tstr},
             missing=["",'""',"NA"],
             impute=True,
             key='eid',
             force = True
             )
    mt = mt.annotate_cols(pheno=ht[mt.s])
    mt = mt.annotate_rows(pheno_mac=hl.struct())

    for phenotype in phenotypes:
        
        # deep copy
        tmt = mt
        tmt = tmt.filter_cols(hl.is_defined(tmt.pheno[phenotype]))
        #tmt = tmt.annotate_rows(MAC = hl.min(hl.agg.call_stats(tmt.GT, tmt.alleles).AC))
        tmt = tmt.annotate_rows(MAC = hl.agg.sum(tmt.DS))

        # annotate original MatrixTable
        mt = mt.transmute_rows(pheno_mac = mt.pheno_mac.annotate(placeholder = tmt.rows()[mt.row_key].MAC))
        mt = mt.transmute_rows(pheno_mac = mt.pheno_mac.rename({"placeholder" : phenotype}))

    # write MAC by phenotype
    ht = mt.select_rows(mt.rsid, mt.pheno_mac).rows()
    ht.flatten().export(out_prefix + ".txt.gz")

    # rename to info field 
    mt = mt.rename({"pheno_mac" : "info" })
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='')
    parser.add_argument('--input_type', default=None, help='')
    parser.add_argument('--out_prefix', default=None, help='')
    parser.add_argument('--out_type', default=None, help='')
    parser.add_argument('--pheno_path', default=None, help='')
    parser.add_argument('--phenotypes', default=None, action=SplitArgsByComma, help='String of phenotypes to be processed')
    
    args = parser.parse_args()

    main(args)

