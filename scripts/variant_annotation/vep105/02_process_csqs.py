#!/usr/bin/env python3

import hail as hl
import argparse
import pandas
import os

from gnomad.utils.vep import process_consequences
from ukb_utils import hail_init
from ko_utils import io

def convert_revel_scores_to_unique_float(table, field='revel_score'):
    # Define a transformation function
    def transform_string_to_float_array(s):
        parts = hl.str(s).split(',')

        float_parts = parts.map(
            lambda x: hl.if_else((x != ".") & hl.is_defined(
                hl.float64(x)), hl.float64(x), hl.null('float64'))
        ).filter(lambda x: hl.is_defined(x))

        return float_parts.first()

    # Apply the transformation function to the revel_score field
    new_table = table.annotate(
        vep=table.vep.annotate(
            transcript_consequences=table.vep.transcript_consequences.map(
                lambda tc: tc.annotate(
                    revel_score=transform_string_to_float_array(tc[field])
                )
            )
        )
    )

    return new_table


def main(args):

    # parser
    input_path = args.input_path
    out_prefix = args.out_prefix

    # standard revel usage
    ht = hl.read_table(input_path)
    ht = convert_revel_scores_to_unique_float(ht, 'revel_score_float') # we use modification to dbNSFP to extract this
    ht = process_consequences(ht)
    ht.write(out_prefix + ".ht", overwrite=True)
    ht.export(out_prefix + ".txt.gz")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--out_prefix', default=None,help='Path prefix for output dataset')

    args = parser.parse_args()

    main(args)

