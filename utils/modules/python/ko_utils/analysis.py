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

def annotation_case_builder(worst_csq_expr, use_loftee: bool = True):
    r'''Annotate consequence categories for downstream analysis
    :param worst_csq_by_gene_canonical_expr: A struct that should contain "most_severe_consequence"
    :param use_loftee: if True will annotate PTVs as either high confidence (ptv) or low confidence (ptv_LC)

    '''
    case = hl.case(missing_false=True)
    if use_loftee:
        case = (case
                .when(worst_csq_expr.lof == 'HC', 'ptv')
                .when(worst_csq_expr.lof == 'LC', 'ptv_LC')
                )
    else:
        case = case.when(
            hl.set(PLOF_CSQS).contains(
                worst_csq_by_gene_canonical_expr.most_severe_consequence),
            "ptv")
    case = (case.when(hl.set(MISSENSE_CSQS).contains(worst_csq_expr.most_severe_consequence) &
                      (~hl.is_defined(worst_csq_expr.cadd_phred) | ~hl.is_defined(worst_csq_expr.revel_score)), "other_missense")
                .when(hl.set(MISSENSE_CSQS).contains(worst_csq_expr.most_severe_consequence) &
                      (worst_csq_expr.cadd_phred >= 20) & (worst_csq_expr.revel_score >= 0.6), "damaging_missense")
                .when(hl.set(MISSENSE_CSQS).contains(worst_csq_expr.most_severe_consequence), "other_missense")
                .when(hl.set(SYNONYMOUS_CSQS).contains(worst_csq_expr.most_severe_consequence), "synonymous")
                .when(hl.set(OTHER_CSQS).contains(worst_csq_expr.most_severe_consequence), "non_coding")
            )
    return case.or_missing()

def count_urv_by_samples(mt):
    r'''Count up ultra rare variants by cols and cateogry
    :param mt: a MatrixTable with the field "consequence_category"
    '''
    return mt.annotate_cols(n_coding_URV_SNP = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_snp(mt.alleles[0], mt.alleles[1]) & (mt.consequence_category != "non_coding")),
                          n_coding_URV_indel = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_indel(mt.alleles[0], mt.alleles[1]) & (mt.consequence_category != "non_coding")),
                          n_URV_PTV = hl.agg.count_where(mt.GT.is_non_ref() & (mt.consequence_category == "ptv")),
                          n_URV_PTV_LC = hl.agg.count_where(mt.GT.is_non_ref() & (mt.consequence_category == "ptv_lc")),
                          n_URV_damaging_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.consequence_category == "damaging_missense")),
                          n_URV_other_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.consequence_category == "other_missense")),
                          n_URV_synonymous = hl.agg.count_where(mt.GT.is_non_ref() & (mt.consequence_category == "synonymous")),
                          n_URV_non_coding = hl.agg.count_where(mt.GT.is_non_ref() & (mt.consequence_category == "non_coding"))
                         )

def count_urv_by_genes(mt, call_gene_expr):
    r''' Count up URVs by gene (as defined by vep.worst_csq_for_variant_canonical.gene_id) '''
    return (mt.group_rows_by(call_gene_expr).
            aggregate(gene = hl.struct(
                n_coding_URV_SNP = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_snp(mt.alleles[0], mt.alleles[1]) & (mt.consequence_category != "non_coding")),
                n_coding_URV_indel = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_indel(mt.alleles[0], mt.alleles[1]) & (mt.consequence_category != "non_coding")),
                n_URV_PTV = hl.agg.count_where(mt.GT.is_non_ref() & (mt.consequence_category == "ptv")),
                n_URV_PTV_LC = hl.agg.count_where(mt.GT.is_non_ref() & (mt.consequence_category == "ptv_LC")),
                n_URV_damaging_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.consequence_category == "damaging_missense")),
                n_URV_other_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.consequence_category == "other_missense")),
                n_URV_synonymous = hl.agg.count_where(mt.GT.is_non_ref() & (mt.consequence_category == "synonymous")),
                n_URV_non_coding = hl.agg.count_where(mt.GT.is_non_ref() & (mt.consequence_category == "non_coding"))
            )      
        )
    )

def annotate_phased_entries(mt):
    r'''Annotates alleles that have the alternate allele on either first or second strand.'''
    mt = mt.annotate_entries(a0_alt=mt.GT == hl.parse_call('1|0'))
    mt = mt.annotate_entries(a1_alt=mt.GT == hl.parse_call('0|1'))
    mt = mt.annotate_entries(a_homo=mt.GT == hl.parse_call('1|1'))
    return mt

def gene_csqs_knockout_builder(in_mt, keep=None):
    '''Return a hail table that contains knockout status alongside phase of variants in genes'''
    mt_rs = gene_strand_builder(in_mt, 'snpid')
    mt_dt = gene_csqs_case_builder(in_mt)
    combined = mt_rs.annotate(csqs=mt_dt[mt_rs.gene_id, mt_rs.s].csqs)
    if keep is not None:
        combined = combined.filter(hl.literal(keep).contains(combined.csqs))
    return combined


def gene_burden_annotations_per_sample(mt):
    r'''count non-ref genotypes per sample'''
    mt = mt.group_rows_by(
        gene_id=mt.consequence.vep.worst_csq_by_gene_canonical.gene_id
    ).aggregate(n=hl.agg.count_where(mt.GT.is_non_ref()))
    return mt


def gene_burden_category_annotations_per_sample(mt):
    r'''count non-ref genotypes per consequence category per sample'''
    mt = mt.group_rows_by(
        gene_id=mt.consequence.vep.worst_csq_by_gene_canonical.gene_id,
        consequence_category=mt.vep.consequence_category
    ).aggregate(n=hl.agg.count_where(mt.GT.is_non_ref()))
    return mt


def gene_burden_stats_per_sample(mt, gt_stats_expr):
    r'''count summary per gene per sample'''
    mt = mt.group_rows_by(
        gene_id=mt.consequence.vep.worst_csq_by_gene_canonical.gene_id,
        consequence_category=mt.vep.consequence_category
    ).aggregate(n=hl.agg.count_where(gt_stats_expr))
    return mt


def maf_category_case_builder(call_stats_expr):
    return (hl.case()
            .when(call_stats_expr.AF <= 0.00001, 0.00001)
            .when(call_stats_expr.AF <= 0.0001, 0.0001)
            .when(call_stats_expr.AF <= 0.001, 0.001)
            .when(call_stats_expr.AF <= 0.01, 0.01)
            .when(call_stats_expr.AF <= 0.1, 0.1)
            .default(0.99))


def mac_category_case_builder(call_stats_expr):
    return (hl.case()
            .when(call_stats_expr.AC <= 5, call_stats_expr.AC)
            .when(call_stats_expr.AC <= 10, 10)
            .when(call_stats_expr.AC <= 25, 25)
            .when(call_stats_expr.AC <= 100, 100)
            .when(call_stats_expr.AC <= 1000, 1000)
            .when(call_stats_expr.AC <= 10000, 10000)
            .when(call_stats_expr.AC <= 100000, 100000)
            .when(call_stats_expr.AC <= 1000000, 1000000)
            .default(0))

def create_gene_map_ht(ht, check_gene_contigs=False):
    #from gnomad.utils.vep import process_consequences
    #ht = process_consequences(ht)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    ht = ht.annotate(
        variant_id=ht.locus.contig + ':' + hl.str(ht.locus.position) + '_' + ht.alleles[0] + '/' + ht.alleles[1],
        annotation=annotation_case_builder(ht.vep.worst_csq_by_gene_canonical)
        )
    if check_gene_contigs:
        gene_contigs = ht.group_by(
            gene_id=ht.vep.worst_csq_by_gene_canonical.gene_id,
            gene_symbol=ht.vep.worst_csq_by_gene_canonical.gene_symbol,
        ).aggregate(
            contigs=hl.agg.collect_as_set(ht.locus.contig)
        )
        assert gene_contigs.all(hl.len(gene_contigs.contigs) == 1)

    gene_map_ht = ht.group_by(
        gene_id=ht.vep.worst_csq_by_gene_canonical.gene_id,
        gene_symbol=ht.vep.worst_csq_by_gene_canonical.gene_symbol,
    ).partition_hint(100).aggregate(
        interval=hl.interval(
            start=hl.locus(hl.agg.take(ht.locus.contig, 1)[0], hl.agg.min(ht.locus.position)),
            end=hl.locus(hl.agg.take(ht.locus.contig, 1)[0], hl.agg.max(ht.locus.position))
        ),
        variants=hl.agg.group_by(ht.annotation, hl.agg.collect(ht.variant_id)),
    )
    return gene_map_ht


def post_process_gene_map_ht(gene_ht):
    groups = ['pLoF', 'damaging_missense|LC', 'pLoF|damaging_missense|LC', 'pLoF|damaging_missense',  'damaging_missense', 'other_missense', 'synonymous']
    variant_groups = hl.map(lambda group: group.split('\\|').flatmap(lambda csq: gene_ht.variants.get(csq)), groups)
    gene_ht = gene_ht.transmute(
        variant_groups=hl.zip(groups, variant_groups)
    ).explode('variant_groups')
    gene_ht = gene_ht.transmute(annotation=gene_ht.variant_groups[0], variants=hl.sorted(gene_ht.variant_groups[1]))
    gene_ht = gene_ht.key_by(start=gene_ht.interval.start)
    return gene_ht.filter(hl.len(gene_ht.variants) > 0)

def count_variants(vep_ht_path, vep_vcf_path):
    from gnomad.utils.vep import process_consequences
    ht = hl.read_table(vep_ht_path)
    ht = process_consequences(ht)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    ht = ht.annotate(
        variant_id=ht.locus.contig + ':' + hl.str(ht.locus.position) + '_' + ht.alleles[0] + '/' + ht.alleles[1],
        annotation=annotation_case_builder(ht.vep.worst_csq_by_gene_canonical)
        )
    ht = ht.filter(hl.literal({'ptv', 'ptv_LC', 'damaging_missense', 'other_missense', 'synonymous'}).contains(ht.annotation))
    print(ht.count())



