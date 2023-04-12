#!/usr/bin/env python3

import hail as hl
import argparse
import random
import string

from ukb_utils import hail_init
from ukb_utils import tables
from ko_utils import io

# Rename the alleles in mt2_shifted using 4 random characters
#def random_string(length):
#    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=length))


def main(args):
    
    chrom = args.chrom
    recessive_path = args.recessive_path
    recessive_type = args.recessive_type
    additive_path = args.additive_path
    additive_type = args.additive_type
    out_prefix = args.out_prefix
    out_type = args.out_type
   
    hail_init.hail_bmrc_init('logs/hail/combine_recessive_common.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1')
 
    rec = io.import_table(recessive_path, recessive_type, calc_info=False)
    add = io.import_table(additive_path, additive_type, calc_info=False)

    # ensure same type before mapping
    add = add.transmute_entries(DS = hl.int32(add.DS))

    # Create a new row with variants by combining chromosome, position, and alleles using ":" as a delimiter
    add = add.annotate_rows(
        rsid = hl.delimit(
            [hl.str(add.locus.contig),
             hl.str(add.locus.position),
             add.alleles[0],
             add.alleles[1]],
            ':')
        ) 

    #print(add_shifted.describe())
    #print(add_shifted.show())
    #print(rec.describe())
    #print(rec.show())
    print(add.rsid.show())

    # Write a table of the new row with variants
    new_variants_table = add.select_rows(add.rsid).rows()
    new_variants_table.export(out_prefix + '.tsv', delimiter='\t')

    # Join the two MatrixTables
    combined_mt = rec.union_rows(add)
    io.export_table(combined_mt, out_prefix, out_type)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='')
    parser.add_argument('--additive_path', default=None, help='')
    parser.add_argument('--additive_type', default=None, help='')
    parser.add_argument('--recessive_path', default=None, help='')
    parser.add_argument('--recessive_type', default=None, help='')
    parser.add_argument('--out_prefix', default=None, help='')
    parser.add_argument('--out_type', default=None, help='')
    args = parser.parse_args()

    main(args)

