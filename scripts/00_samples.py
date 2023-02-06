#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import samples

def main(args):

    output_path = args.output_path
    export_imputed_samples = args.export_imputed_samples
    export_genotyped_samples = args.export_genotyped_samples 
    get_unrelated = args.get_unrelated
    remove_withdrawn = args.remove_withdrawn
    extract1 = args.extract1
    extract2 = args.extract2

    if export_imputed_samples:
        mt = genotypes.get_ukb_imputed_v3_bgen([21])

    if export_genotyped_samples:
        mt = genotypes.get_ukb_genotypes_bed([21])

    if remove_withdrawn:
        mt = samples.remove_withdrawn(mt)

    if extract1:
        ht_samples = hl.import_table(
            extract1,
            no_header=True, key='f0',
            delimiter=',')
        mt = mt.filter_cols(
            hl.is_defined(ht_samples[mt.col_key]))

    if extract2:
        ht_samples = hl.import_table(
            extract2,
            no_header=True, key='f0',
            delimiter=',')
        mt = mt.filter_cols(
            hl.is_defined(ht_samples[mt.col_key]))

    # filter to unrelated
    if get_unrelated:
        mt = samples.filter_ukb_to_unrelated_using_kinship(mt)

    # export
    mt.s.export(output_path, header = False)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_path', default=None, help='')
    parser.add_argument('--export_imputed_samples', action='store_true', default=None, help='')
    parser.add_argument('--export_genotyped_samples', action='store_true', default=None, help='')
    parser.add_argument('--get_unrelated', action='store_true', default=None, help='')
    parser.add_argument('--remove_withdrawn', action='store_true', default=None, help='')
    parser.add_argument('--extract1', default=None, help='')
    parser.add_argument('--extract2', default=None, help='')

    args = parser.parse_args()

    main(args)



