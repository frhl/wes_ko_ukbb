#!/usr/bin/env python3

from ko_utils import readphasing as rp

import argparse
import numpy as np
import pandas as pd


def main(args):

    path_fasta_reference = args.path_fasta_reference
    path_informative_snps = args.path_informative_snps
    path_query_snps = args.path_query_snps
    chromosome = args.chromosome
    outfile = args.outfile

    # load all potential informative variants (slow)
    df_full = pd.read_csv(path_informative_snps, sep = '\t')
        
    # subset to chromosome in question
    df_query = pd.read_csv(path_query_snps, sep = '\t')
    df_query = df_query.loc[df_query.chromosome == chromosome]
    df_query = df_query.reset_index()  # make sure indexes pair with number of rows
    # get list of list of target variants
    target_variants = np.array(rp.query_target_variants(df_query, exclude_phased = True))
    query_rows = np.shape(target_variants)[0]

    # ensure that target variants are in df_full
    reference_variants = df_full.varid.tolist()
    variants_in_reference = np.array([var in reference_variants for var in target_variants[:,1].tolist()])
    assert np.all(variants_in_reference) 

    print("Iterating over %s rows for chromosome %s using %s" % (
        query_rows, chromosome, outfile
        ))
    with open(outfile,"w") as outfile:
        for idx in range(query_rows):

            # extract target
            row = target_variants[idx,:]
            eid = int(row[0])
            target_variant = row[1]
            target_genotype = row[2]

            # get paths and positions
            target_position = int(target_variant.split(":")[1])
            cram_path = rp.get_cram_path(eid)

            # bit of verbose information
            print("%s with %s (%s)" % (eid, target_variant, target_genotype))
            print("Running with CRAM file: %s" % (cram_path))

            # get informative varinats
            informative_variants = rp.get_informative_variants(df_full, eid, target_variant, 5)
            r = rp.phase(eid, 
                      chromosome, 
                      target_position, 
                      target_variant, 
                      informative_variants, 
                      cram_path, 
                      path_fasta_reference
                     )
            
            outstring = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                r["eid"],
                r["target_variant"],
                r["target_phase_prediction"],
                r["query_variant"],
                r["query_phase"],
                r["distance"],
                r["n_overlaps"],
                r["class_prediction"],
                r["evidence_count"],
                r["cis_evidence"],
                r["trans_evidence"]
            )
            outfile.write(outstring)
        

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--path_fasta_reference', default=None, help='fasta file with reference genome')
    parser.add_argument('--path_informative_snps', default=None, help='path to file containing all informative SNPs, including SNPs that should be queried')
    parser.add_argument('--path_query_snps', default=None, help='path to file of SNPs to query')
    parser.add_argument('--outfile', default=None, help='fasta file with reference genome')
    parser.add_argument('--chromosome', default=None, help='fasta file with reference genome')
    
    args = parser.parse_args()

    main(args)

