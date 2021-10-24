#!/usr/bin/env python3

def flip_base(base: str) -> str:
    r""" 
    :param base: a string that is is in ATGC
    """    
    return (hl.switch(base)
                .when('A','T')
                .when('T','A')
                .when('G','C')
                .when('C','G')
                .default(base))


def liftover(mt: MatrixTable, current = 'GRCh37', to = 'GRCh38', drop_annotations = False) -> MatrixTable:
    r""" Liftover variants from one build to another

    :param mt: a Matrixtable
    :param current: the current build (either 'GRCh37' or 'GRCh38')
    :param to: the desired build
    :param drop_annotations: boolean indicating whether old locus and alleles should be dropped. 
    """
    
    builds = ['GRCh37','GRCh38']
    assert(current in builds)
    assert(to in builds)    

    # setup liftover references
    liftover_path = variants.get_liftover_chain_path(current,to)
    rg_from = hl.get_reference(current)  
    rg_to = hl.get_reference(to)  
    rg_from.add_liftover(liftover_path, rg_to)

    # do liftover and flip allele for negative strand
    mt = mt.annotate_rows(new_locus=hl.liftover(mt.locus, to, include_strand=True),old_locus=mt.locus)  
    mt = mt.filter_rows(hl.is_defined(mt.new_locus))
    mt = mt.annotate_rows(new_alleles=hl.if_else(mt.new_locus.is_negative_strand, [flip_base(mt.alleles[0]), flip_base(mt.alleles[1])],mt.alleles))
    mt = mt.key_rows_by(locus=mt.new_locus.result, alleles=mt.new_alleles)

    # drop old fields
    if drop_annotations:
        mt = mt.drop('new_locus')
        mt = mt.drop('old_locus')
        mt = mt.drop('new_alleles')

    return mt
 

