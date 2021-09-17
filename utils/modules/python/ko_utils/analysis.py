#!/usr/bin/env python3

import hail as hl

def annotate_vep(mt, vep_path):
    r'''Annotate matrix table with VEP consequence from external file.'''
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
    
    # Most severe variant consequence
    ht = ht.annotate_rows(vep = ht.vep.annotate(most_severe_consequence = ht.vep.Consequence.split('&')[0]))
    
    # Extract various categories annotations and change type
    ht = ht.annotate_rows(vep = ht.vep.annotate(sift_pred = ht.vep.SIFT_pred.split('&')[0]))
    ht = ht.annotate_rows(vep = ht.vep.annotate(polyphen2_hdiv_pred = ht.vep.Polyphen2_HDIV_pred.split('&')[0]))
    ht = ht.annotate_rows(vep = ht.vep.annotate(polyphen2_hvar_pred = ht.vep.Polyphen2_HVAR_pred.split('&')[0]))
    ht = ht.annotate_rows(vep = ht.vep.annotate(cadd_phred_score = hl.parse_float(ht.vep.CADD_phred)))
    ht = ht.annotate_rows(vep = ht.vep.annotate(revel_score = hl.parse_float(ht.vep.REVEL_score)))
    
    # Define protein truncating variants
    ptv = hl.set(["transcript_ablation", "splice_acceptor_variant",
              "splice_donor_variant", "stop_gained", "frameshift_variant"])
    
    # Define missense variation
    missense = hl.set(["stop_lost", "start_lost", "transcript_amplification",
                   "inframe_insertion", "inframe_deletion", "missense_variant",
                   "protein_altering_variant", "splice_region_variant"])
    
    # Define synonymous
    synonymous = hl.set(["incomplete_terminal_codon_variant", "stop_retained_variant", "synonymous_variant"])
    
    # Define non coding variation
    non_coding = hl.set(["coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant",
              "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant"])
    
    # Create categories for downstream analysis
    ht = ht.annotate_rows(vep = ht.vep.annotate(consequence_category = 
        hl.case().when(ptv.contains(ht.vep.most_severe_consequence), "ptv")
             .when(missense.contains(ht.vep.most_severe_consequence) & 
                   (~hl.is_defined(ht.vep.cadd_phred_score) | 
                    ~hl.is_defined(ht.vep.revel_score)), "other_missense")                                   
             .when(missense.contains(ht.vep.most_severe_consequence) & 
                   (ht.vep.cadd_phred_score >= 20) & 
                   (ht.vep.revel_score >= 0.6), "damaging_missense") 
             .when(missense.contains(ht.vep.most_severe_consequence), "other_missense")
             .when(synonymous.contains(ht.vep.most_severe_consequence), "synonymous")
             .when(non_coding.contains(ht.vep.most_severe_consequence), "non_coding")
             .default("NA")
    ))
                                                
    # combine with matrix table
    mt = mt.annotate_rows(vep = ht.index_rows(mt.locus, mt.alleles).vep.drop('CSQ'))
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
    # get all snps that are not homozygous
    mt = in_mt
    mt = analysis.annotate_phased_entries(mt)
    mt = mt.filter_entries(~mt.GT.is_hom_var())
    # create table for each strand and combine to gene
    ht0 = (mt.group_rows_by(mt.vep.Gene).aggregate(s0 = hl.agg.any(mt.a0_alt)))
    ht1 = (mt.group_rows_by(mt.vep.Gene).aggregate(s1 = hl.agg.any(mt.a1_alt)))
    ht2 = (in_mt.group_rows_by(in_mt.vep.Gene).aggregate(hom_var = hl.agg.any(in_mt.GT.is_hom_var())))
    # combine entries
    ht = ht0.annotate_entries(s1 = ht1[ht0.Gene, ht0.s].s1)
    ht = ht.annotate_entries(hom_var = ht2[ht.Gene, ht.s].hom_var)
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
    # get all snps that are not homozygous
    mt = in_mt
    mt = analysis.annotate_phased_entries(mt)
    mt = mt.filter_entries(~mt.GT.is_hom_var())
    # create table for each strand and combine to gene
    ht0 = (mt.group_rows_by(mt.vep.Gene).aggregate(s0 = hl.agg.any(mt.a0_alt)))
    ht1 = (mt.group_rows_by(mt.vep.Gene).aggregate(s1 = hl.agg.any(mt.a1_alt)))
    ht2 = (in_mt.group_rows_by(in_mt.vep.Gene).aggregate(hom_var = hl.agg.any(in_mt.GT.is_hom_var())))
    # combine entries
    ht = ht0.annotate_entries(s1 = ht1[ht0.Gene, ht0.s].s1)
    ht = ht.annotate_entries(hom_var = ht2[ht.Gene, ht.s].hom_var)
    expr = (hl.case()
           .when( (ht.s0) & (ht.s1) & (ht.hom_var), 2)
           .when( (ht.s0) & (ht.s1), 2)
           .when( (ht.hom_var), 2)
           .when( (ht.s0) & (ht.s1 == False), 1)
           .when( (ht.s1) & (ht.s0 == False), 1)
           .default(0))
    ht = ht.annotate_entries(DT = expr)
    ht = ht.drop('s0').drop('s1').drop('hom_var')    
    return ht


def count_alleles(mt):
    r'''Count up alleles in a vector (singleton AC, not singletons AC, total AC)'''
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
    d = mt.aggregate_entries(hl.agg.group_by(mt.vep.Gene, hl.agg.count_where(mt.GT.is_non_ref())))
    n_genes = len([(key, value) for key, value in d.items() if value != 0])
    return n_genes

def summarize_variants(mt, what = 'ptv', vep_field = 'consequence_category'):
    r'''Count up singletons, non singletons, total and genes affected by variants'''
    ht = mt.filter_rows(hl.literal(set([what])).contains(mt.vep[vep_field]))
    ht_alleles = count_alleles(ht)
    ht_genes = count_genes(ht)
    out = ht_alleles + [ht_genes]
    return out

def gene_burden_annotations_per_sample(mt, gene_field = 'Gene'):
    r''' calculate gene burden by counting variants in gene'''
    mt = mt.group_rows_by(
        Gene = mt.vep.Gene#,
        #consequence_category = mt.vep.consequence_category
        ).aggregate(n = hl.agg.count_where(mt.GT.is_non_ref()))
    return mt

def gene_burden_category_annotations_per_sample(mt, gene_field = 'Gene'):
    r''' calculate gene burden by counting variants in gene stratified by variant category'''
    mt = mt.group_rows_by(
        Gene = mt.vep.Gene,
        consequence_category = mt.vep.consequence_category
        ).aggregate(n = hl.agg.count_where(mt.GT.is_non_ref()))
    return mt

def calc_p_ko(mt):
    '''Annotates entries with P(Knockout). Requires, that fields are
       already annotated with "singletons" count, and "DT". '''
    ko_mt = mt.annotate_entries(
        pKO = hl.if_else(
            mt.DT == 2, 1, # knockout
            hl.if_else(
                mt.DT == 1, 
                hl.if_else(mt.singletons >= 1, 1 - (1/2)**mt.singletons, 0), # one phased hetz
                hl.if_else(mt.singletons >= 2, 1 - 2*(1/2)**mt.singletons, 0), # zero phased hetz
            )
        )
    )
    return ko_mt

def gene_csqs_calc_pKO(mt_phased, mt_unphased, fields_drop = ['dosage','sigletons']):
    
    '''Annotates entries with P(Knockout). Requires a phased matrix and an unphased matrix that only contains singletons.'''
    
    # setup variables
    mt1 = mt_phased # contains phased non-singletons
    mt2 = mt_unphased # contains unphased singletons
    
    # Determine probability of being KO given singletons and phased hetz
    mt1_burden = gene_csqs_dosage_builder(mt1)
    mt2_burden = gene_burden_annotations_per_sample(mt2)
    mt_ko = mt1_burden.annotate_entries(singletons = mt2_burden[(mt1_burden.Gene, mt1_burden.s)].n)
    mt_ko = mt_ko.annotate_entries(singletons = hl.if_else(~hl.is_missing(mt_ko.singletons), mt_ko.singletons, 0 ))
    mt_ko = calc_p_ko(mt_ko)

    # drop not needed rows
    mt_ko = mt_ko.drop(*[f for f in fields_drop if f in mt_ko.entry])
    return mt_ko

def gene_csqs_calc_pKO_pseudoSNP(mt1, mt2, chrom):
    '''Calculate probability of being a knockout incoporating phased 
       data (mt1) and unphased singletons (mt2). Create a file with 
       fake markers, that can be inputted into SAIGE'''
    # mt1 is phased
    # mt2 is unphased

    # get probability matrix
    pmt = get_prob_ko_matrix(mt1, mt2, ["DT","singletons"])
    pmt = pmt.annotate_entries(DT = pmt.pKO*2) # multiply probability by 2 (dosage encoded [0:2])
    pmt = pmt.drop('pKO')

    # create fake loci
    pmt = pmt.annotate_rows(locus = hl.parse_locus('chr' + str(chrom) + ':1'))
    pmt = pmt.annotate_rows(alleles = hl.literal(['X','Y']))
    pmt = pmt.annotate_rows(rsid = pmt.Gene)
    pmt = pmt.key_rows_by(pmt.locus, pmt.alleles)
    pmt = pmt.drop('Gene')
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

