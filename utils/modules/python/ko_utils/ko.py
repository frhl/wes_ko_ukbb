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


def csqs_case_builder(worst_csq_expr: hl.StringExpression, use_loftee: bool = True, loftee_lc_annotation="damaging_missense"):
    r'''Annotate consequence categories for downstream analysis
    
    :param worst_csq_by_gene_canonical_expr: A struct that should contain "most_severe_consequence"
    :param use_loftee: if True will annotate PTVs as either high confidence (ptv) or low confidence (ptv_LC)
    :param loftee_lc_annotation: low confidence
    '''
    case = hl.case(missing_false=True)
    if use_loftee:
        assert loftee_lc_annotation in ("LC", "damaging_missense")
        print(f"Note: LOFTEE Low Confidence PTVs annotated as {loftee_lc_annotation}")
        case = (case
                .when(worst_csq_expr.lof == 'HC', 'pLoF')
                .when(worst_csq_expr.lof == 'LC', loftee_lc_annotation)
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


def unphase(gt: hl.call):
    """ unphase a genotype call 
    
    :param gt: GT (call)
    """
    assert str(gt.dtype) == 'call'
    return((hl.case()
    .when(gt == hl.parse_call("1|0"), hl.parse_call("0/1"))
    .when(gt == hl.parse_call("0|1"), hl.parse_call("0/1"))
    .when(gt == hl.parse_call("1|1"), hl.parse_call("1/1"))
    .when(gt == hl.parse_call("0|0"), hl.parse_call("0/0"))
    .when(gt == hl.parse_call("1/0"), hl.parse_call("0/1"))
    .when(gt == hl.parse_call("0/1"), hl.parse_call("0/1"))
    .when(gt == hl.parse_call("1/1"), hl.parse_call("1/1"))
    .when(gt == hl.parse_call("0/0"), hl.parse_call("0/0"))
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


def rand_hom_to_het(gt: hl.call, P: float = 0.9, seed = None, phase = "|"):
    """ Randomize genotype phase of call
    :param gt: genotype call to be flipped
    :param P: probabily of flipping hom to het
    :param seed: seed for random
    """
    assert str(gt.dtype) == 'call'

    return hl.if_else(
        (gt.n_alt_alleles() == 2),
        hl.if_else(
            (hl.rand_bool(P, seed=seed)),
            hl.if_else(
                (hl.rand_bool(0.5, seed=seed)),
                hl.parse_call("1" + str(phase) + "0"),
                hl.parse_call("0" + str(phase) + "1")
            ),
            gt
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

   

def calc_prob_ko(hom_expr, phased_expr, unphased_expr, only_homs=False):
    """Calculate probability of knockout based on phased 
       and unphased alleles (requires sum_gts_entries)
    
    :param hom_expr: integer for homozygous count 
    :param phased_expr: struct with integers a1 and a2
    :param unphased_expr: struct with integers a1 and a2
    :param only_homs: only count homozygotes as knockouts.
    """
    
    n = unphased_expr.n
    #n = unphased_expr.a1 + unphased_expr.a2 
    if only_homs:
        return (hl.case()
               .when(hom_expr > 0, 1) # homozygote
               .default(0)
               )
    else:
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


def calc_prob_ko_by_count(ko_expr, phased_expr, unphased_expr):
    """Calculate probability of knockout based on a knockout expression
    and a count of phased and unphased het sites. Unlike calc_prop_ko,
    this function is designed to take in an expression containing counts
    and from this calculate the probability of the locus being knocked out.
    
    :param ko_expr: integer for knockout count (either homozygous
    or compound heterozygous knockouts)
    :param phased_expr: count of phased heterozygous sites
    :param unphased_expr: count of unphased heterozygous sites
    """
    
    return (hl.case()
        .when(
            ko_expr == 1, 1)
        .when(
            (phased_expr == 1) & (unphased_expr > 0), 
            (1 - (1 / 2) ** unphased_expr))
        .when(
            (phased_expr == 0) & (unphased_expr > 1), 
            (1 - 2 * (1 / 2) ** unphased_expr))
        .default(0))



def annotate_knockout(hom_expr, pko_expr, phased_expr = None):
    """Annotate entry knockout type 
   
    :param hom_expr: integer for homozygous count 
    :param pko_expr: probability of knockout (see calc_prob_ko)
    :param phased_expr: struct with integers a1 and a2. 
    """
    if phased_expr:
        return (hl.case()
               .when((hom_expr > 0) & (pko_expr == 1), 'Homozygote')
               .when((hom_expr == 0) & (pko_expr == 1),'Compound heterozygote')
               .when((hom_expr == 0) & (pko_expr >= 0.5),'Possible Compound heterozygote')
               .when((hom_expr == 0) & (pko_expr == 0) & (phased_expr.a1 > 1),'Compound heterozygote (cis)')
               .when((hom_expr == 0) & (pko_expr == 0) & (phased_expr.a2 > 1),'Compound heterozygote (cis)')
               .when((hom_expr == 0) & (pko_expr == 0) & (phased_expr.a1 == 1),'Heterozygote')
               .when((hom_expr == 0) & (pko_expr == 0) & (phased_expr.a2 == 1),'Heterozygote')
               .or_missing()
                )
    else:
        return (hl.case()
               .when((hom_expr > 0) & (pko_expr == 1), 'Homozygote')
               .when((hom_expr == 0) & (pko_expr == 1),'Compound heterozygote')
               .when((hom_expr == 0) & (pko_expr >= 0.5),'Possible Compound heterozygote')
               .or_missing()
                )

def discard_prob_dosages(DS):
    """remove any compound het owed to unphased
    singletons. Will return either 0 or 2.
    :param ds: DOSAGE float between 0 and 2. 
    """
    return (hl.case().when(DS == 2, 2).default(0))



def normalize_by_name(mt, name):
    """Normalize entry to have mean of zero and variance of 1"""
    mt = mt.annotate_rows(**{'stats': hl.agg.stats(mt[name])})
    mt = mt.annotate_entries(
        norm=(mt[name]-mt.stats.mean)/mt.stats.stdev)
    mt = mt.rename({"norm" : "norm_" + str(name)})
    return(mt)


def make_thetas(mt, h2, pi = None):
    """ Make thetas (effect sizes) for genes (gene x sample matrix) with
    either infintesimal or spike and slab model. 
    """
    
    M = mt.count()[0]
    pi_temp = 1 if pi == None else pi
    mt = mt.annotate_rows(
        theta = hl.rand_bool(pi_temp)*hl.rand_norm(0, hl.sqrt(h2/(M*pi_temp))))
    return(mt)

def simulate_effect_size(mt, h2, pi = None):
    """ Simulate effect sizes for either variants/genes under 
    the infintesimal or spike and slab model. 
    :param mt: MatrixTable
    :param h2: Heritability
    :param pi: probability of a gene/variant being causal
    """
    
    M = mt.count()[0]
    pi_temp = 1 if pi == None else pi
    return(hl.rand_bool(pi_temp)*hl.rand_norm(0, hl.sqrt(h2/(M*pi_temp))))


def make_effect_size(mt, beta, pi = None):
    """ Make specified effect sizes for genes (gene x sample matrix)
    regardless of final heritability.
    :param mt: MatrixTable
    :param beta: float for pre-specified effect size.   
    :param pi: probability that a gene is causal.
    """
    
    pi_temp = 1 if pi == None else pi
    return(hl.rand_bool(pi_temp)*beta)

def get_gt_from_floor_ds(DS):
    """ Take dosage and create fake GT calls. This ensures that VCFs
    with GT annotations can be read back in via hail. Note that DS is
    rounded down to nearest integer (0, 1, 2).
    :param DS: Dosage (float64)
    """
    return (hl.case()
     .when(hl.int(hl.floor(DS)) == 0, hl.parse_call("0/0"))
     .when(hl.int(hl.floor(DS)) == 1, hl.parse_call("1/0"))
     .when(hl.int(hl.floor(DS)) == 2, hl.parse_call("1/1"))
     .default(hl.parse_call("0/0")))
 




