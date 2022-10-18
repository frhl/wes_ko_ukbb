#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init

def chunks(lst, n):
    """Yield successive n-sized chunks from lst.
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def write_interval_file(samples, path, n):
    """Write interval file to path
    :param samples: list of IDs of samples (e.g. from mt.s.collect())
    :param path: file to write to
    :param n: how big should each chunk be
    """
    assert n > 0
    with open(path, "w") as outfile:
        for chunk in chunks(samples, n):
            outline = ",".join(chunk) + "\n"
            outfile.write(outline)


def read_interval_idx(interval_path, idx):
    """ Read a specific interval index from a file
    :param interval_path: path to interval file (assumine one line interval)
    :param idx: interval index
    """
    assert idx >= 0
    with open(interval_path, "r") as infile:
        for i, line in enumerate(infile):
            if i == idx:
                return(line)


def main(args):
    
    input_path = args.input_path
    output_path = args.output_path
    input_type = args.input_type
    output_type = args.output_type
    samples_per_chunk = args.samples_per_chunk
    interval_path = args.interval_path
    interval_idx = args.interval_idx
    write_interval = args.write_interval
    split_by_interval = args.split_by_interval
    

    hail_init.hail_bmrc_init_local('logs/hail/prephase_chunks.log', 'GRCh38')
    mt = io.import_table(input_path, input_type, calc_info = False)
    
    if write_interval:
        samples = mt.s.collect()
        write_interval_file(samples, interval_path, int(samples_per_chunk))
    elif split_by_interval:
        if interval_idx:
            samples = read_interval_idx(interval_path, int(interval_idx) - 1)
            samples = samples.strip().split(",")
            mt = mt.filter_cols(hl.literal(samples).contains(mt.s))
            io.export_table(mt, output_path, output_type)
        else:
            raise TypeError("'split_by_interval' specified but not 'phasing_idx'!")
    else:
        raise ValueError("No command given")




if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', default=None, help='What is the input directory of the files?')
    parser.add_argument('--output_path', default=None, help='What is the input directory of the files?')
    parser.add_argument('--input_type', default=".vcf.gz", help='What extension does the file(s) end with?')
    parser.add_argument('--output_type', default=".vcf.gz", help='What extension does the file(s) end with?')
    parser.add_argument('--samples_per_chunk', default=None, help='What extension does the file(s) end with?')
    parser.add_argument('--interval_path', default=None, help='How big should the new overlap be?')
    parser.add_argument('--interval_idx', default=None, help='current interval number in 1-based index')
    parser.add_argument('--write_interval', default=None, action='store_true', help='should interval file be written')
    parser.add_argument('--split_by_interval', default=None, action='store_true', help='should mt be split by interval_idx')

    args = parser.parse_args()

    main(args)


