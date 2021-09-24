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

def variant_csqs_category_builder(mt):
    r'''Create categories for downstream analysis'''
    return mt.annotate_rows(vep = mt.vep.annotate(consequence_category = 
        hl.case().when(hl.literal(set(PLOF_CSQS)).contains(mt.vep.most_severe_consequence), "ptv")
             .when(hl.literal(set(MISSENSE_CSQS)).contains(mt.vep.most_severe_consequence) & 
                   (~hl.is_defined(mt.dbnsfp.cadd_phred_score) | 
                    ~hl.is_defined(mt.dbnsfp.revel_score)), "other_missense")                                   
             .when(hl.literal(set(MISSENSE_CSQS)).contains(mt.vep.most_severe_consequence) & 
                   (mt.dbnsfp.cadd_phred_score >= 20) & 
                   (mt.dbnsfp.revel_score >= 0.6), "damaging_missense") 
             .when(hl.literal(set(MISSENSE_CSQS)).contains(mt.vep.most_severe_consequence), "other_missense")
             .when(hl.literal(set(SYNONYMOUS_CSQS)).contains(mt.vep.most_severe_consequence), "synonymous")
             .when(hl.literal(set(OTHER_CSQS)).contains(mt.vep.most_severe_consequence), "non_coding")
             .default("NA")))

def annotate_dbnsfp(mt, vep_path):
    r'''Annotate matrix table with dbNSFP consequence from external VEP file.'''
    print(f'Annotating with VEP file: {vep_path}')
    
    # Open file containing VEP fields
    with open('data/vep/vep_fields.txt', 'r') as file:
        fields = file.read().strip().split(',')
    ht = hl.import_vcf(vep_path).rename({'info':'vep'}) 
    
    # Add VEP fields by iteration
    for i in range(len(fields)):
        ht = ht.annotate_rows(
            vep=ht.vep.annotate(
                col=ht.vep.CSQ.map(lambda x: (x.split('\\|')[i]))[0]
                ).rename({'col':f'{fields[i]}'})
        )
    
    # Extract various categories annotations and change type
    ht = ht.annotate_rows(vep = ht.vep.annotate(sift_pred = ht.vep.SIFT_pred.split('&')[0]))
    ht = ht.annotate_rows(vep = ht.vep.annotate(polyphen2_hdiv_pred = ht.vep.Polyphen2_HDIV_pred.split('&')[0]))
    ht = ht.annotate_rows(vep = ht.vep.annotate(polyphen2_hvar_pred = ht.vep.Polyphen2_HVAR_pred.split('&')[0]))
    ht = ht.annotate_rows(vep = ht.vep.annotate(cadd_phred_score = hl.parse_float(ht.vep.CADD_phred)))
    ht = ht.annotate_rows(vep = ht.vep.annotate(revel_score = hl.parse_float(ht.vep.REVEL_score)))
    
    # annotate main table
    mt = mt.annotate_rows(dbnsfp = hl.struct())
    mt = mt.annotate_rows(dbnsfp = mt.dbnsfp.annotate(revel_score = ht.index_rows(mt.locus, mt.alleles).vep.revel_score))
    mt = mt.annotate_rows(dbnsfp = mt.dbnsfp.annotate(cadd_phred_score = ht.index_rows(mt.locus, mt.alleles).vep.cadd_phred_score))
    mt = mt.annotate_rows(dbnsfp = mt.dbnsfp.annotate(polyphen2_hdiv_pred = ht.index_rows(mt.locus, mt.alleles).vep.polyphen2_hdiv_pred))
    mt = mt.annotate_rows(dbnsfp = mt.dbnsfp.annotate(polyphen2_hvar_pred = ht.index_rows(mt.locus, mt.alleles).vep.polyphen2_hvar_pred))

    return(mt)



def filter_vep(mt, field, conds):
    r'''Filter VEP field by condition(s) '''
    mt = mt.filter_rows(hl.literal(set(conds)).contains(mt.vep[field]))
    return mt


def annotate_phased_entries(mt):
    r'''Annotates alleles that have the alternate allele on either first or second strand.'''
    mt = mt.annotate_entries(a0_alt = mt.GT ==  hl.parse_call('1|0'))
    mt = mt.annotate_entries(a1_alt = mt.GT ==  hl.parse_call('0|1'))
    mt = mt.annotate_entries(a_homo = mt.GT ==  hl.parse_call('1|1'))
    return mt

def gene_csqs_case_builder(in_mt):
    r''' Returns matrix table that contains gene consequence information from phased geneotypes.
    "": no alternate alleles,
    "HE": one alternate allele on either strand in a locus, 
    "HO": homozygous for alternate alleles
    "CH": two alternate allele on either strand in a locus (compound heterozygous)
    "CH+HO": two alternate allele on either strand in a locus (either as homozygous or compound heterozygous)
    '''
    # create one gene_id for each item in gene_id array
    in_mt = in_mt.explode_rows(in_mt.vep.worst_csq_by_gene_canonical)
    # get all snps that are not homozygous
    mt = in_mt
    mt = annotate_phased_entries(mt)
    mt = mt.filter_entries(~mt.GT.is_hom_var())
    # create table for each strand and combine to gene
    ht0 = (mt.group_rows_by(mt.vep.worst_csq_by_gene_canonical.gene_id).aggregate(s0 = hl.agg.any(mt.a0_alt)))
    ht1 = (mt.group_rows_by(mt.vep.worst_csq_by_gene_canonical.gene_id).aggregate(s1 = hl.agg.any(mt.a1_alt)))
    ht2 = (in_mt.group_rows_by(in_mt.vep.worst_csq_by_gene_canonical.gene_id).aggregate(hom_var = hl.agg.any(in_mt.GT.is_hom_var())))
    # combine entries
    ht = ht0.annotate_entries(s1 = ht1[ht0.gene_id, ht0.s].s1)
    ht = ht.annotate_entries(hom_var = ht2[ht.gene_id, ht.s].hom_var)
    expr = (hl.case()
           .when( (ht.s0) & (ht.s1) & (ht.hom_var), 'CH+HO')
           .when( (ht.s0) & (ht.s1), "CH")
           .when( (ht.hom_var), 'HO')
           .when( (ht.s0) & (ht.s1 == False), 'HE')
           .when( (ht.s1) & (ht.s0 == False), 'HE')
           .default(''))
    ht = ht.annotate_entries(csqs = expr)
    ht = ht.drop('s0').drop('s1').drop('hom_var')    
    return ht


def gene_csqs_dosage_builder(in_mt):
    r''' Returns matrix table that contains dosage information from phased geneotypes.
    0: two refererence alleles in locus,
    1: one alternate allele on either strand in a locus, 
    2: two alternate allele on either strand in a locus (either as homozygous or compound heterozygous)
    '''
    # create one gene_id for each item in gene_id array
    in_mt = in_mt.explode_rows(in_mt.vep.worst_csq_by_gene_canonical)
    # get all snps that are not homozygous
    mt = in_mt
    mt = annotate_phased_entries(mt)
    mt = mt.filter_entries(~mt.GT.is_hom_var())
    # create table for each strand and combine to gene
    ht0 = (mt.group_rows_by(mt.vep.worst_csq_by_gene_canonical.gene_id).aggregate(s0 = hl.agg.any(mt.a0_alt)))
    ht1 = (mt.group_rows_by(mt.vep.worst_csq_by_gene_canonical.gene_id).aggregate(s1 = hl.agg.any(mt.a1_alt)))
    ht2 = (in_mt.group_rows_by(in_mt.vep.worst_csq_by_gene_canonical.gene_id).aggregate(hom_var = hl.agg.any(in_mt.GT.is_hom_var())))
    # combine entries
    ht = ht0.annotate_entries(s1 = ht1[ht0.gene_id, ht0.s].s1)
    ht = ht.annotate_entries(hom_var = ht2[ht.gene_id, ht.s].hom_var)
    expr = (hl.case()
           .when( (ht.s0) & (ht.s1) & (ht.hom_var), 2)
           .when( (ht.s0) & (ht.s1), 2)
           .when( (ht.hom_var), 2)
           .when( (ht.s0) & (ht.s1 == False), 1)
           .when( (ht.s1) & (ht.s0 == False), 1)
           .default(0))
    ht = ht.annotate_entries(DS = expr)
    ht = ht.drop('s0').drop('s1').drop('hom_var')    
    return ht

def count_alleles(mt):
    r'''Count up alleles in a vector (singleton AC, not singletons AC, total AC)'''
    mt = mt.explode_rows(mt.vep.worst_csq_by_gene_canonical)
    mt = mt.filter_rows(mt.info.AC > 0)
    d = mt.aggregate_rows(hl.agg.counter(mt.info.AC))
    if 1 in d.keys():
        n_singletons_ac = d[1]
    else:
        n_singletons_ac = 0
    n_not_singletons_ac = sum([value*key for key, value in d.items() if key != 1])
    return [n_singletons_ac, n_not_singletons_ac, n_singletons_ac + n_not_singletons_ac]

def count_genes(mt):
    r'''Collapse variants into genes and count affected genes'''
    mt = mt.explode_rows(mt.vep.worst_csq_by_gene_canonical)
    d = mt.aggregate_entries(hl.agg.group_by(mt.vep.worst_csq_by_gene_canonical.gene_id, hl.agg.count_where(mt.GT.is_non_ref())))
    n_genes = len([(key, value) for key, value in d.items() if value != 0])
    return n_genes

def summarize_variants(mt, what = 'ptv', vep_field = 'consequence_category'):
    r'''Count up singletons, non singletons, total and genes affected by variants'''
    ht = mt.filter_rows(hl.literal(set([what])).contains(mt.vep[vep_field]))
    ht_alleles = count_alleles(ht)
    ht_genes = count_genes(ht)
    out = ht_alleles + [ht_genes]
    return out

def gene_burden_annotations_per_sample(mt):
    r''' calculate gene burden by counting variants in gene'''
    mt = mt.explode_rows(mt.vep.worst_csq_by_gene_canonical)
    mt = mt.group_rows_by(
        gene_id = mt.vep.worst_csq_by_gene_canonical.gene_id
        ).aggregate(n = hl.agg.count_where(mt.GT.is_non_ref()))
    return mt

def gene_burden_category_annotations_per_sample(mt):
    r''' calculate gene burden by counting variants in gene stratified by variant category'''
    mt = mt.explode_rows(mt.vep.worst_csq_by_gene_canonical)
    mt = mt.group_rows_by(
        gene_id = mt.vep.worst_csq_by_gene_canonical.gene_id,
        consequence_category = mt.vep.consequence_category
        ).aggregate(n = hl.agg.count_where(mt.GT.is_non_ref()))
    return mt

def calc_p_ko(mt):
    '''Annotates entries with P(Knockout). Requires, that fields are
       already annotated with "singletons" count, and "DS". '''
    ko_mt = mt.annotate_entries(
        pKO = hl.if_else(
            mt.DS == 2, 1, # knockout
            hl.if_else(
                mt.DS == 1, 
                hl.if_else(mt.singletons >= 1, 1 - (1/2)**mt.singletons, 0), # one phased hetz
                hl.if_else(mt.singletons >= 2, 1 - 2*(1/2)**mt.singletons, 0), # zero phased hetz
            )
        )
    )
    return ko_mt

def gene_csqs_calc_pKO(mt_phased, mt_unphased, fields_drop = ['DS','singletons']):
    '''Annotates entries with P(Knockout). Requires a phased matrix and an unphased matrix that only contains singletons.'''
    
    # setup variables
    mt1 = mt_phased # contains phased non-singletons
    mt2 = mt_unphased # contains unphased singletons
    
    # Determine probability of being KO given singletons and phased hetz
    mt1_burden = gene_csqs_dosage_builder(mt1)
    mt2_burden = gene_burden_annotations_per_sample(mt2)
    mt_ko = mt1_burden.annotate_entries(singletons = mt2_burden[(mt1_burden.gene_id, mt1_burden.s)].n)
    mt_ko = mt_ko.annotate_entries(singletons = hl.if_else(~hl.is_missing(mt_ko.singletons), mt_ko.singletons, 0 ))
    mt_ko = calc_p_ko(mt_ko)

    # drop not needed rows
    if fields_drop is not None:
        mt_ko = mt_ko.drop(*[f for f in fields_drop if f in mt_ko.entry])
    return mt_ko


def gene_csqs_calc_pKO_pseudoSNP(mt1, mt2, chrom):
    '''Calculate probability of being a knockout incoporating phased 
       data (mt1) and unphased singletons (mt2). Create a file with 
       fake markers, that can be inputted into SAIGE'''
    # mt1 is phased
    # mt2 is unphased

    # get probability matrix
    pmt = gene_csqs_calc_pKO(mt1, mt2, ["DS","singletons"])
    pmt = pmt.annotate_entries(DS = pmt.pKO*2) # multiply probability by 2 (dosage encoded [0:2])
    pmt = pmt.drop('pKO')

    # create fake loci
    pmt = pmt.annotate_rows(locus = hl.parse_locus('chr' + str(chrom) + ':1'))
    pmt = pmt.annotate_rows(alleles = hl.literal(['X','Y']))
    pmt = pmt.annotate_rows(rsid = pmt.gene_id)
    pmt = pmt.key_rows_by(pmt.locus, pmt.alleles)
    pmt = pmt.drop('gene_id')
    return pmt

def maf_category_case_builder(mt):
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

def gene_strand_builder(mt, field = 'snpid'):
    '''Returns hail table that contains genes, samples, rsids, knockout status'''
    
    # annotate entries with phased data
    mt = mt.explode_rows(mt.vep.worst_csq_by_gene_canonical)
    mt = mt.annotate_entries(rsid_entry = mt[field])
    mt = mt.annotate_rows(gene_id = mt.vep.worst_csq_by_gene_canonical.gene_id)
    mt = annotate_phased_entries(mt)

    # filter to each strand
    strand1 = mt.filter_entries((mt.a0_alt == True) | (mt.a_homo == True)).entries() # sets entries to NA in matrix table
    strand2 = mt.filter_entries((mt.a1_alt == True) | (mt.a_homo == True)).entries()

    # filter to each gene
    strand1 = strand1.group_by(strand1.gene_id, strand1.s).aggregate(phase1 = hl.agg.collect(strand1.rsid_entry))
    strand2 = strand2.group_by(strand2.gene_id, strand2.s).aggregate(phase2 = hl.agg.collect(strand2.rsid_entry))

    # combine each strand
    ht = strand1.annotate(phase2 = strand2[strand1.gene_id, strand1.s].phase2)
    ht = ht.annotate(knockout = (~hl.is_missing(ht.phase1)) & ~(hl.is_missing(ht.phase2)))
    return ht

def gene_csqs_knockout_builder(in_mt, keep = None):
    '''Return a hail table that contains knockout status alongside phase of variants in genes'''
    mt_rs = gene_strand_builder(in_mt, 'snpid')
    mt_dt = gene_csqs_case_builder(in_mt)
    combined = mt_rs.annotate(csqs = mt_dt[mt_rs.gene_id, mt_rs.s].csqs)
    if keep is not None:
        combined = combined.filter(hl.literal(keep).contains(combined.csqs))
    return combined

