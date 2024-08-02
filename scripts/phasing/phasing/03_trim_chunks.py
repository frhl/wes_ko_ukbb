#!/usr/bin/env python3

import numpy as np
import hail as hl
import argparse
import math
import os
import re

from ko_utils import io
from ukb_utils import hail_init


def sort_by_chunks(chunks):
    ''' sort by X in XofY with regex '''
    n = len(chunks)
    n_chunks_expt = int(re.findall(r'\d+of\d+', chunks[0])[0].split("of")[1])
    assert n == n_chunks_expt
    lst = { int(re.findall(r'\d+of\d+', c)[0].split("of")[0])-1 : c for c in chunks }
    chunks = [lst[i] for i in range(n)]
    return chunks

def main(args):
    
    in_dir = args.in_dir
    in_ext = args.in_ext
    new_overlap_size = args.new_overlap_size
    in_prefix = args.in_prefix
    out_prefix = args.out_prefix
    out_type = args.out_type

    # set up
    hail_init.hail_bmrc_init_local('logs/hail/merge_chunks.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip

    # assume directory contains only phased vcf files
    files = [os.path.join(in_dir, f) for f in os.listdir(in_dir) if f.endswith(in_ext)]
    if len(files) == 0:
        raise ValueError(f"No files with extension {in_ext} found at {in_dir}")
    if in_prefix:
        files = [f for f in files if in_prefix in f]
    if len(files) == 0:
        raise ValueError("Files were found. But no file with input prefix were found. Exiting!")
    
    # import each table as a matrix-table:
    print(files)
    files = sort_by_chunks(files)
    mts = [hl.import_vcf(f, force_bgz = True, find_replace=(':-nan', ':NaN')) for f in files]
    n = len(mts)

    with open(out_prefix + "_trims.txt", "w") as outfile:
        if len(mts) > 1:
            i = 0
            max_iter = len(mts)-1
            mt1 = mts[0]
            mt2 = mts[1]
            final = list()
            while (len(final) < max_iter):
            
                # get overlapping fraction
                overlap = mt1.filter_rows(hl.is_defined(mt2.rows()[mt1.row_key]))
                n_overlap = overlap.count()[0]

                # desired left and right flank 
                overlap = overlap.add_row_index()
                right_flank_end = math.floor(n_overlap / 2) + math.floor(int(new_overlap_size) / 2)
                left_flank_start = math.ceil(n_overlap / 2) - math.ceil(int(new_overlap_size) / 2)

                # append indexes to overlap, and select right flank for mt1 and left flank for mt2
                mt1_right_flank = overlap.filter_rows(
                        overlap.row_idx == right_flank_end ).locus.position.collect()[0]
                mt2_left_flank = overlap.filter_rows(
                        overlap.row_idx == left_flank_start).locus.position.collect()[0]
                
                # get position for start/end of overlap
                start = overlap.head(1).locus.position.collect()[0]
                end = overlap.tail(1).locus.position.collect()[0]
                outfile.write(f"{in_dir}\t{in_prefix}\t{i}\t{mt1_right_flank}\t{mt2_left_flank}\t{overlap}\t{start}\t{end}\n")
                
                # Trim right flank.
                mt1 = mt1.filter_rows(mt1.locus.position <= mt1_right_flank)
                
                # if mt1 is the first matrix (from left to right) it is now done.
                # otherwise continue and cut down left flank.
                final.append(mt1)
                
                # Finish left flank.
                # if mt2 is the last matrix, it's now done.
                # if mt2 is the middle matrix, the left flank is now done. 
                # Set mt2 to mt1 in order to trim its right flank
                mt2 = mt2.filter_rows(mt2.locus.position >= mt2_left_flank)
                
                i += 1
                mt1 = mt2
                if (i < max_iter):
                    mt2 = mts[i+1]
            
            # add right-most matrix table
            final.append(mt2)
            
            # write out the trimmed chunks
            j = 0
            for mt in final:
                j += 1
                out = out_prefix + "_trim." + str(j) + 'of' + str(n) 
                outfile = out + ".vcf.bgz"
                if not os.path.isfile(outfile):
                    print(f"Writing {out}.vcf.bgz ..")
                    io.export_table(mt, out, "vcf")

        else:
            mt = mts[0]
            out = out_prefix + "_trim.1of1" 
            outfile = out + ".vcf.bgz"
            if not os.path.isfile(outfile):
                print(f"Writing {out}.vcf.bgz ..")
                io.export_table(mt, out, out_type)

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--in_dir', default=None, help='What is the input directory of the files?')
    parser.add_argument('--in_ext', default=".vcf.gz", help='What extension does the file(s) end with?')
    parser.add_argument('--in_prefix', default=None, help='What extension does the file(s) end with?')
    parser.add_argument('--new_overlap_size', default=None, help='How big should the new overlap be?')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default="vcf", help='vcf/plink/mt')
    args = parser.parse_args()

    main(args)


