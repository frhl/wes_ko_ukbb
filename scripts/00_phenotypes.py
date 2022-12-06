#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init

def main(args):

    # parser
    input_path = args.input_path
    extract_samples = args.extract_samples
    extract_phased_samples = args.extract_phased_samples
    export_header = args.export_header
    count_case_control = args.count_case_control
    only_males = args.only_males
    only_females = args.only_females
    out_prefix = args.out_prefix

    # setup flags
    hail_init.hail_bmrc_init_local('logs/hail/phenotypes.log', 'GRCh38')

    ht = hl.import_table(input_path, impute = True, key = 'eid', missing = ["NA",""," "], types = {"eid": hl.tstr}, force = True)

    if extract_samples:
        samples = hl.import_table(extract_samples,no_header=True, key='f0',delimiter=',')
        ht = ht.filter(hl.is_defined(samples[ht.key]))

    if extract_phased_samples:
        phased_samples = hl.import_table(extract_phased_samples, no_header=True, key='f0')
        ht = ht.filter(hl.is_defined(phased_samples[ht.key]))

    if only_females:
        ht = ht.filter(ht.sex == 1)
    
    if only_males:
        ht = ht.filter(ht.sex == 0)

    ht.export(out_prefix + ".tsv")

    # quick and dirty subset to
    # the columns that only contain phenotypes
    cols = list(ht.row)
    cols = cols[10:len(cols)]
    phenos = [c for c in cols if "PC" not in c]
    ht = ht.select(*phenos)

    if export_header:
        outfile = out_prefix + "_header.tsv"
        with open(outfile, "w") as outfile:
            for pheno in phenos:
                outfile.write(pheno + "\n")

    if count_case_control:
        outfile = out_prefix + "_counts.tsv"
        print(ht.describe())
        with open(outfile, "w") as outfile:
            for pheno in phenos:
                cases = ht.aggregate(hl.agg.sum(ht[pheno] == True))
                controls = ht.aggregate(hl.agg.sum(ht[pheno] == False))
                # avoid division by zero
                if (cases + controls) > 0:
                    fraction = cases / (cases + controls)
                else:
                    fraction = 0
                line = ("%s\t%d\t%d\t%f" % (pheno, cases, controls, fraction))  
                outfile.write(line + "\n")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_path', default=None, help='Path to input phenotypes')
    parser.add_argument('--extract_samples', default=None, help='Path to dbNSFP annotations')
    parser.add_argument('--extract_phased_samples', default=None, help='Path to dbNSFP annotations')
    parser.add_argument('--export_header', default=None, action='store_true', help='Path to dbNSFP annotations')
    parser.add_argument('--count_case_control', default=None, action='store_true', help='Path to dbNSFP annotations') 
    parser.add_argument('--only_males', default=None, action='store_true',  help='Path to dbNSFP annotations')
    parser.add_argument('--only_females', default=None, action='store_true', help='Path to dbNSFP annotations')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')

    args = parser.parse_args()

    main(args)



