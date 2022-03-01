#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init

def main(args):

    # parser
    input_path = args.input_path
    extract_samples = args.extract_samples
    only_males = args.only_males
    only_females = args.only_females
    out_prefix = args.out_prefix

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip

    ht = hl.import_table(input_path, impute = True, key = 'eid', missing = ["NA",""], types = {"eid": hl.tstr})

    if extract_samples:
        samples = hl.import_table(extract_samples,no_header=True, key='f0',delimiter=',')
        ht = ht.filter(hl.is_defined(samples[ht.key]))

    if only_females:
        ht = ht.filter(ht.sex == 1)
    
    if only_males:
        ht = ht.filter(ht.sex == 0)

    ht.export(out_prefix + ".tsv.gz")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_path', default=None, help='Path to input phenotypes')
    parser.add_argument('--extract_samples', default=None, help='Path to dbNSFP annotations')
    parser.add_argument('--only_males', default=None, help='Path to dbNSFP annotations')
    parser.add_argument('--only_females', default=None, help='Path to dbNSFP annotations')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')

    args = parser.parse_args()

    main(args)



