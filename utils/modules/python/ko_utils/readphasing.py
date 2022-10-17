#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import pysam

def alelles_at_read_position(read, position):
    """ Get allele at current select position in the read.
    :param read: the current read from samfile.fetch
    :param position: position based on reference coordinates
    """
    reference_positions = read.get_reference_positions()
    sequence = list(read.query_alignment_sequence)
    position_index = reference_positions.index(position-1)
    return(sequence[position_index])

def get_reads_supported_by_allele_snv(chrom, position, cram, reference):
    """ Get dictionary of reads that support the allele at the position
    :param chrom: chromosome input
    :param pos: position in the reference genome
    :param cram: cram file containing aligned reads
    :param reference: reference genome fasta file
    """
    handle = pysam.AlignmentFile(cram, "rc", reference_filename = reference)
    reads_by_allele = {base : [] for base in "ACGTN"}
    for read in handle.fetch(chrom, pos-1, pos):
        allele = alelles_at_read_position(read, position)
        reads_by_allele[allele].append(read.query_name)
    handle.close()
    return(reads_by_allele)

def get_alleles_supported_by_read(chrom, position, cram, reference):
    """ Get dictionary of alleles that are supported by reads
    :param chrom: chromosome input
    :param pos: position in the reference genome
    :param cram: cram file containing aligned reads
    :param reference: reference genome fasta file
    """
    handle = pysam.AlignmentFile(cram, "rc", reference_filename = reference)
    allele_by_reads = {}
    for read in handle.fetch(chrom, position-1, position):
        try:
            allele = alelles_at_read_position(read, position)
            allele_by_reads[read.query_name] = allele 
        except ValueError:
            allele_by_reads[read.query_name] = None
    handle.close()
    return(allele_by_reads)

def get_informative_variants(df, eid, variant, n_nearby = 3, exclude_target = True, exclude_unphased = True):
    """ Return N variants nearest to target variant based on reference chromosome position. This
    assume that the df has already been ordered by position, which is default from Hail.
    :param df: a data.frame containing all samples/snps
    :param eid: a sample ID that is also in df
    :param variant: string in the format of chr:pos:a1:a2
    :param n_nearby: how many snps upstream/downstream should be retrieved?
    :param exclude_target: should target snp be excluded from the resulting list?
    :parma exclude unphased: should unphased variants be excluded? This will also remove the target variant
    """

    # deep copy to avoid changing main dataframe
    df_subset = df.copy()
    df_subset = df_subset.loc[df_subset.s == int(eid)]
    n_rows = df_subset.shape[0]
    # where is the singleton located?
    idx = np.where(df_subset.varid == variant)[0].tolist()[0]
    # get n nearest heterozygous alleles upstream and downstream
    idx_lower_bound = max(0, idx-n_nearby)
    idx_upper_bound = min(n_rows, idx+n_nearby+1)
    # get variants around target variant
    target_variant = df_subset.iloc[idx].varid
    nearby_variants = df_subset.iloc[idx_lower_bound:idx_upper_bound]
    # get variants and resulting genotypes
    variants = nearby_variants.varid.tolist()
    genotypes = nearby_variants.GT.tolist()
    # exclude target variant?
    if exclude_target:
        idx_to_pop = variants.index(target_variant)
        variants.pop(idx_to_pop)
        genotypes.pop(idx_to_pop)
    if exclude_unphased:
        idx_keep = [ i for (i,x) in enumerate(genotypes) if "/" not in x]
        variants = [variants[i] for i in idx_keep]
        genotypes = [genotypes[i] for i in idx_keep]
    # create dict of variant to genotype
    assert len(genotypes) > 1
    out = {variants[i] : genotypes[i] for i in range(len(variants))}
    return(out)
    
def get_is_alt(variant):
    """ Build dict to check if inputted allele is alternate with respect to the variant
    """
    splitted = variant.split(":")
    return(
        dict({ 
            splitted[2] : False,
            splitted[3] : True
        })
    )

def subset_df_to_singletons(df):
    df_singletons = df.copy()
    idx = np.where(
        (df_singletons.GT == "0/1") | 
        (df_singletons.GT == "1/0")
    )[0].tolist()
    return(df_singletons.iloc[idx])


def phase(eid, chrom, target_position, target_variant, variants_to_query, cram_path, ref):
    """ Perform readbacked phasing between a target snp and a query variants
    
    """

    # get current target allele
    target_support = get_alleles_supported_by_read(chrom, target_position, cram_path, ref)
    target_reads_id = target_support.keys()
    target_alleles_id = target_support.values()
    target_is_alt = get_is_alt(target_variant)

    for var in enumerate(variants_to_query.keys()):    
        
        # get query reads
        query_variant = str(var[1])
        query_position = int(var[1].split(":")[1])
        query_support = get_alleles_supported_by_read(chrom, query_position, cram_path, ref)
        query_reads_id = query_support.keys()
        query_alleles = query_support.values()
        query_is_alt = get_is_alt(query_variant)
        query_phase = variants_to_query[query_variant]
        distance = abs(query_position-target_position)

        # calculate overlapping reads
        read_overlap = list(set(target_reads_id) & set(query_reads_id))
        n_overlaps = len(read_overlap)
        read_backed_evidence = {"cis" : 0, "trans" : 0}
    
        # init variables
        target_phase_prediction = None
        class_prediction = None
        evidence_count = 0 
            
        
        if n_overlaps > 0:

            for idx in range(n_overlaps):

                # what is the allele at the target position?
                target_allele_in_read = target_support[read_overlap[idx]]
                target_allele_is_alt = target_is_alt[target_allele_in_read]

                # what is the allele at the informative/query position?
                query_allele_in_read = query_support[read_overlap[idx]]
                query_allele_is_alt = query_is_alt[query_allele_in_read]


                # check whether cis or trans evidence by read
                if target_allele_is_alt and query_allele_is_alt:
                    read_backed_evidence["cis"] += 1
                elif target_allele_is_alt and not query_allele_is_alt:
                    read_backed_evidence["trans"] += 1
                elif not target_allele_is_alt and query_allele_is_alt:
                    read_backed_evidence["trans"] += 1
                elif not target_allele_is_alt and not query_allele_is_alt:
                    read_backed_evidence["cis"] += 1
                else:
                    raise TypeError("Unexpected combination of REF/ALT")

        # make a prediction based on most votes from cis or trans
        if read_backed_evidence["cis"] == read_backed_evidence["trans"]:
            class_prediction = "unknown"
        elif read_backed_evidence["cis"] >= read_backed_evidence["trans"]:
            evidence_count = read_backed_evidence["cis"]
            target_phase_prediction = query_phase
            class_prediction = "cis"

        else:
            evidence_count = read_backed_evidence["trans"]
            target_phase_prediction = query_phase[::-1]
            class_prediction = "trans"

        # print to stdout or file
        out = {
            "eid" : eid,
            "target_variant" : target_variant,
            "target_phase_prediction" : target_phase_prediction,
            "query_variant" : query_variant,
            "query_phase" : query_phase,
            "distance" : distance,
            "n_overlaps" : n_overlaps,
            "class_prediction" : class_prediction,
            "evidence_count" : evidence_count,
            "cis_evidence" : read_backed_evidence["cis"],
            "trans_evidence" : read_backed_evidence["trans"]
            }

        return(out)

    


def query_target_variants(df, exclude_phased = False):
    """ Query variants from data.frame and return in list of lists containing
    sample ids, varaints, and the corresponding ohase
    :param df: data.frame with column s, varid (csv) and gts (csv)
    :param exclude_phase: boolean indicate whether phased hits should be excluded
    """

    assert "varid" in df
    assert "gts" in df
    assert "s" in df
    lst = list()
    for index, row in df_query.iterrows():
        eid = row['s']
        variants = row['varid'].split(",")
        phases = row['gts'].split(",")
        for idx in range(len(variants)):
            current_variant = variants[idx]
            current_phase = phases[idx]
            append_to_list = exclude_phased and "|" in current_phase
            if not append_to_list:
                lst.append([
                    eid, 
                    current_variant,
                    current_phase,
                ])
    return(lst)

def get_cram_path(eid):
    """Return CRAM file for a specific sample """
    path = "/well/ukbb-wes/cram/oqfe/ukbb-11867/" + str(eid) + "_oqfe.cram"
    if os.path.exists(path):
        return(path)
    else:
        raise TypeError(str(path), + " does not exists!")
    


