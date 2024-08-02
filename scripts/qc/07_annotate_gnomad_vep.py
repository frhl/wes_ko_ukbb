#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import qc
from ko_utils import analysis
from gnomad.utils.vep import process_consequences

def main(args):
    
    # parser
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    #mt = qc.get_table(input_path=input_path, input_type=input_type)
    mt = hl.import_vcf(input_path, force_bgz = True)

    # subset to rows
    mt = hl.vep(mt, 'utils/configs/vep_env.json')
    mt = process_consequences(mt)
    ht = mt.rows().select('vep')
    
    ht.write(out_prefix + "_vep_qc.ht", overwrite=True)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')


    args = parser.parse_args()

    main(args)

