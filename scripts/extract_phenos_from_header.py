#!/usr/bin/env python3

# A script that extract a list of strings from a header (e.g. a header).

from ko_utils import phenos
import argparse
import sys

def main(args):
    
    out = phenos.extract_from_header(args.input, delim = args.input_delim)
    if args.index is not None:
        out = out[int(args.index)-1]
        sys.stdout.write(out)
    else:
        sys.stdout.write("\t".join(out))

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default=None, help='Path to phenotype file')
    parser.add_argument('--regex', default=None, help='What regex should be used?')
    parser.add_argument('--index', default=None, help='select an index 0-Inf to be carried forward')
    parser.add_argument('--input_delim', default="\t", help='delimiter for input file')
    args = parser.parse_args()
    main(args)

