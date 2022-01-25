#!/usr/bin/env python3

import hail as hl


def filter_max_af(mt, af=None):
    r'''Filter to variants to have af less than {maf}'''
    mt = mt.filter_rows(hl.min(mt.info.AF) < af)
    return mt


def filter_min_af(mt, af=None):
    r'''Filter to variants to have af gt {maf}'''
    mt = mt.filter_rows(hl.min(mt.info.AF) > af)
    return mt


def filter_max_maf(mt, maf=None):
    r'''boolean for variants that have maf less than {maf}'''
    maf_expr = (mt.info.AF < maf) | (mt.info.AF > (1 - maf))
    mt = mt.filter_rows(maf_expr)
    return mt


def filter_min_maf(mt, maf=None):
    r'''boolean for variants that have maf gt {maf}'''
    maf_expr = (mt.info.AF > maf) & (mt.info.AF < (1 - maf))
    mt = mt.filter_rows(maf_expr)
    return mt

def filter_maf(mt, min_maf=0, max_maf=0.5):
    r'''Filter to variants based on a certain min/max MAF threshold'''
    assert min_maf >= 0
    assert max_maf <= 0.5
    maf_expr = ((mt.info.AF > min_maf) & (mt.info.AF < (1 - min_maf)))
    maf_expr = maf_expr & ((mt.info.AF < max_maf) | (mt.info.AF > (1 - max_maf)))
    return mt.filter_rows(maf_expr)

def filter_max_mac(mt, mac=None):
    r'''Filter to variants to have maf less than {maf}'''
    if mac is not None:
        mt = mt.filter_rows(hl.min(mt.info.AC) <= mac)
    return mt


def filter_min_mac(mt, mac=None):
    r'''Filter to variants to have mac gt {mac}'''
    if mac is not None:
        mt = mt.filter_rows(hl.min(mt.info.AC) >= mac)
    return mt


def filter_min_missing(mt, value=0.05):
    r'''Filter variants to have less than {value} in genotype missigness"'''
    assert value >= 0 and value < 1
    missing = hl.agg.mean(hl.is_missing(mt.GT)) <= value
    mt = mt.filter_rows(missing)
    return mt


def filter_hwe(mt, cut_off=1e-12):
    r'''Perform filtering based on hwe p-values.'''
    mt = mt.annotate_rows(hwe=hl.agg.hardy_weinberg_test(mt.GT))
    rows = mt.hwe.p_values >= cut_off
    mt = mt.filter_rows(rows)
    return mt


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


