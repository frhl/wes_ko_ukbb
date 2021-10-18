#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import qc

def main(args):
    
    # parser
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    gnomad_path = args.gnomad_path

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    mt = qc.get_table(input_path=input_path, input_type=input_type)

    # subset to rows
    ht = mt.rows()
    ht = hl.vep(ht, vep_config)
    ht = process_consequences(ht)

    # Annotate with gnomAD information for the relevant chromosome
    gnomad_variants_ht = hl.import_vcf(gnomad_path, reference_genome ='GRCh38', force_bgz=True, array_elements_required=False).rows()
    ht = ht.annotate(inGnomAD = hl.is_defined(gnomad_variants_ht[ht.key]))
    #ht = ht.annotate(annotation=annotation_case_builder(ht.vep.worst_csq_for_variant_canonical)) # from ukb_common import *

    ht.write(out_prefix, overwrite=True)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--gnomad_path', default=None, help='Path to gnomAD exomes for comparison"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')


    args = parser.parse_args()

    main(args)

