#!/usr/bin/env python3

import hail as hl

PLOF_CSQS = ["transcript_ablation", "splice_acceptor_variant",
             "splice_donor_variant", "stop_gained", "frameshift_variant"]

MISSENSE_CSQS = ["stop_lost", "start_lost", "transcript_amplification",
                 "inframe_insertion", "inframe_deletion", "missense_variant"]

SYNONYMOUS_CSQS = ["stop_retained_variant", "synonymous_variant"]

OTHER_CSQS = ["mature_miRNA_variant", "5_prime_UTR_variant",
              "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant"]


def csqs_case_builder(worst_csq_expr: hl.StringExpression, use_loftee: bool = True):
    r'''Annotate consequence categories for downstream analysis
    
    :param worst_csq_by_gene_canonical_expr: A struct that should contain "most_severe_consequence"
    :param use_loftee: if True will annotate PTVs as either high confidence (ptv) or low confidence (ptv_LC)
    '''
    case = hl.case(missing_false=True)
    if use_loftee:
        case = (case
                .when(worst_csq_expr.lof == 'HC', 'pLoF')
                .when(worst_csq_expr.lof == 'LC', 'LC')
                )
    else:
        case = case.when(
            hl.set(PLOF_CSQS).contains(
                worst_csq_expr.most_severe_consequence),
            "pLoF")
    case = (case.when(hl.set(MISSENSE_CSQS).contains(worst_csq_expr.most_severe_consequence) &
                      (~hl.is_defined(worst_csq_expr.cadd_phred) | ~hl.is_defined(worst_csq_expr.revel_score)), "other_missense")
                .when(hl.set(MISSENSE_CSQS).contains(worst_csq_expr.most_severe_consequence) &
                      (worst_csq_expr.cadd_phred >= 20) & (worst_csq_expr.revel_score >= 0.6), "damaging_missense")
                .when(hl.set(MISSENSE_CSQS).contains(worst_csq_expr.most_severe_consequence), "other_missense")
                .when(hl.set(SYNONYMOUS_CSQS).contains(worst_csq_expr.most_severe_consequence), "synonymous")
                .when(hl.set(OTHER_CSQS).contains(worst_csq_expr.most_severe_consequence), "non_coding")
            )
    return case.or_missing()


def set_to_phased_call(gt: hl.call):
    """ Set an unphased genotype call to phased genotype 
    
    :param gt: genotype call to be converted
    """
    assert str(gt.dtype) == 'call'
    return((hl.case()
    .when(gt == hl.parse_call("1/0"), hl.parse_call("1|0"))
    .when(gt == hl.parse_call("0/1"), hl.parse_call("0|1"))
    .when(gt == hl.parse_call("1/1"), hl.parse_call("1|1"))
    .when(gt == hl.parse_call("0/0"), hl.parse_call("0|0"))
    .or_missing()))



def is_phased(gt: hl.call):
    """is a call phased
    
    :param gt: GT (call)
    """
    return (
        (gt == hl.parse_call('1|0')) |
        (gt == hl.parse_call('0|1')) |
        (gt == hl.parse_call('1|1')) |
        (gt == hl.parse_call('0|0'))
           )

def rand_flip_call(gt: hl.call, P: float = 0.5, seed = None):
    """ Randomize genotype phase of call

    :param gt: genotype call to be flipped
    :param P: probabily of one phase
    :param seed: seed for random
    """
    assert str(gt.dtype) == 'call'

    return hl.if_else(
        (gt.n_alt_alleles() == 1) &
        (is_phased(gt)), 
        hl.if_else(
            hl.rand_bool(P, seed=seed),
            hl.parse_call("1|0"),
            hl.parse_call("0|1")
        ),
        gt
    )


def aggr_count_calls(mt: hl.MatrixTable, phased: bool = True):
    """ Count number of phased/unphased hetz and what haplotype they reside on

    :param mt: MatrixTable with GT field
    :param mt: test for phased genotypes only
    """
    aggr = mt.aggregate_entries(hl.agg.counter(mt.GT))
    gt10 = aggr[hl.Call(alleles=[1, 0], phased=phased)]
    gt01 = aggr[hl.Call(alleles=[0, 1], phased=phased)]
    return((gt10,gt01))


def collect_phase_count_by_expr(mt: hl.MatrixTable, expr: hl.StringExpression):
    """Create a hail table of aggregated genotypes by expr
    
    :param mt: MatrixTable to be used
    :param expr: what expression to collapse on, e.g. "gene_id"
    """
    return mt.group_rows_by(expr).aggregate(
              gts=hl.agg.filter(mt.GT.is_non_ref(), hl.agg.collect(mt.GT)),
              varid=hl.agg.filter(mt.GT.is_non_ref(), hl.agg.collect(mt.varid))
            )


def aggr_phase_count_by_expr(mt: hl.MatrixTable, expr):
    """Get aggregated phased/unphased. This is a faster method
    than using "sum_gst_entries", however, it does not allow
    for variant IDs to be collected.
    
    :param mt: MatrixTable with at least some phased entries
    :param expr: StringExpression for what string to perform
    the aggregation on.
    
    """
    return (mt.group_rows_by(expr)
                .aggregate(
                    phased=hl.struct(
                        a1=hl.agg.count_where(
                            (mt.GT == hl.parse_call("1|0"))),
                        a2=hl.agg.count_where(
                            (mt.GT == hl.parse_call("0|1"))),
                        n=hl.agg.count_where(
                            (mt.GT.is_het_ref()) & (mt.GT.phased))
                    ),  
                    unphased=hl.struct(
                        n=hl.agg.count_where(
                            (mt.GT.is_het_ref()) & (~mt.GT.phased))
                    ),
                    hom_alt_n=hl.agg.count_where(mt.GT.is_hom_var())
                    )
           )


def sum_gts_entries(mt: hl.MatrixTable):
    """sum phased/unphased calls 

    :param ht: MatrixTable with gts entry (see collect_gt_by_expr)
    """
    assert 'gts' in list(mt.entry), 'missing entry field "gts"'
   
    return (mt.annotate_entries( 
            hom_alt_n=hl.sum(mt.gts.map(
                    lambda x: x.is_hom_var())),
            phased=hl.struct(
                a1=hl.sum(mt.gts.map(
                    lambda x: x == hl.parse_call('1|0'))),
                a2=hl.sum(mt.gts.map(
                    lambda x: x == hl.parse_call('0|1')))
                ),
            unphased=hl.struct(
                n=hl.sum(mt.gts.map(
                    lambda x: (~x.phased) & (x.is_het_ref())
                    )
                )
            )
        )
    )

   

def calc_prob_ko(hom_expr, phased_expr, unphased_expr):
    """Calculate probability of knockout based on phased 
       and unphased alleles (requires sum_gts_entries)
    
    :param hom_expr: integer for homozygous count 
    :param phased_expr: struct with integers a1 and a2
    :param unphased_expr: struct with integers a1 and a2
    """
    
    n = unphased_expr.n
    #n = unphased_expr.a1 + unphased_expr.a2 
    return (hl.case()
           .when(hom_expr > 0, 1) # homozygote
           .when(
               (phased_expr.a1 > 0) & # compound het
               (phased_expr.a2 > 0), 1)
           .when(
               ((phased_expr.a1 == 1) | # likely compound het (one phased het)
                (phased_expr.a2 == 1)) &
                (n > 0), 1 - (1 / 2) ** n)
           .when(
               ((phased_expr.a1 == 0) | # likely compound het (zero phased het)
                (phased_expr.a2 == 0)) &
                (n > 1), 1 - 2 * (1 / 2) ** n)
           .default(0)
            )


def annotate_knockout(hom_expr, pko_expr):
    """Annotate entry knockout type
   
    :param hom_expr: integer for homozygous count 
    :param pko_expr: probability of knockout (see calc_prob_ko)
    """
   
    return (hl.case()
           .when((hom_expr > 0) & (pko_expr == 1), 'Homozygote')
           .when((hom_expr == 0) & (pko_expr == 1),'Compound heterozygote')
           .when((hom_expr == 0) & (pko_expr >= 0.5),'Possible Compound heterozygote')
           .or_missing()
            )


