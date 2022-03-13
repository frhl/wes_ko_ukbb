#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import io
from ko_utils import ko


def main(args):
    
    chrom = args.chrom
    input_path = args.input_path
    input_type = args.input_type
    permutations = args.permutations
    max_permutations = args.max_permutations
    seed = args.seed
    out_prefix = args.out_prefix
    out_type = args.out_type
    
    # run parser
    hail_init.hail_bmrc_init('logs/hail/array_permute.log', 'GRCh38')
 
    # get table of variants involed in knockouts
    ht = hl.import_table(
        permutations,
        force = True,
        impute = True
    )

    # get variants in the speceific chromosome
    ht = ht.filter(ht.CHR == chrom)
    variants = ", ".join(ht.varid.collect()).split(', ') 

    # how many permutations required?
    n_cur = int(ht.aggregate(hl.agg.max(ht.permut)))
    n = min(n_cur, int(max_permutations))
    print(f"{n_cur} permutations required. Using {n}.")

    # get matrix table with variants
    mt = io.import_table(input_path, input_type)
    mt.filter_rows(hl.literal(set(variants)).contains(mt.varid))

    # permute variant phase
    mt_prev = None
    for i in range(n):
        mt_ = mt
        mt_ = mt_.annotate_rows(permute_id=i)
        mt_ = mt_.transmute_entries(
            GT = ko.rand_flip_call(mt_.GT, seed = int(seed)*99)
        )
        if mt_prev is not None:
            mts = mts.union_rows(mt_)
        else:
            mts = mt_
        mt_prev = mt_

    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='')
    parser.add_argument('--input_path', default=None, help='')
    parser.add_argument('--input_type', default=None, help='')
    parser.add_argument('--permutations', default=None, help='')
    parser.add_argument('--max_permutations', default=None, help='')
    parser.add_argument('--seed', default=None, help='Seed used for randomizing')
    parser.add_argument('--out_prefix', default=None, help='')
    parser.add_argument('--out_type', default=None, help='')
    args = parser.parse_args()

    main(args)



