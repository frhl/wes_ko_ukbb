#!/usr/bin/env python3

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



def filter_vep(mt, field, conditions):
    r'''Filter VEP field by condition(s) '''
    assert isinstance(conditions, list)
    assert field in list(mt.vep)
    #conds = [[cond] for cond in conditions]
    mt = mt.filter_rows(hl.literal(set(conds)).contains(mt.vep[field]))
    return mt


def annotate_phased_entries(mt):
    r'''Annotates alleles that have the alternate allele on either first or second strand.'''
    mt = mt.annotate_entries(a0_alt = mt.GT ==  hl.parse_call('1|0'))
    mt = mt.annotate_entries(a1_alt = mt.GT ==  hl.parse_call('0|1'))
    mt = mt.annotate_entries(a_homo = mt.GT ==  hl.parse_call('1|1'))
    return mt


def construct_phased_dosage_mt(mt, gene_field = 'Gene'):
    r''' Returns matrix table that contains dosage information from phased geneotypes.
    0: two refererence alleles in locus,
    1: one alternate allele on either strand in a locus, 
    2: two alternate allele on either strand in a locus (either as homozygous or compound heterozygous)
    '''
    mt = annotate_phased_entries(mt)
    knockout =  (mt.a0_alt & mt.a1_alt) | mt.a_homo # ((mt.a0_alt & mt.a1_alt) | mt.a_homo ) & (mt.vep.consequence_category == hl.literal("ptv"))
    heterozygous = (mt.a0_alt | mt.a1_alt)
    burden_mt = (
        mt 
        .group_rows_by(mt.vep[gene_field])
        .aggregate(dosage = hl.if_else( hl.agg.any(knockout), 2, 
                            hl.if_else( hl.agg.any(heterozygous), 1, 0 )))
        #.burden_mt.filter_entries(burden_mt.entries != 0)
    )
    return burden_mt

def construct_summary_mt(mt, gene_field = 'Gene'):
    r''' Returns matrix table that contains dosage information from phased geneotypes.
    0: two refererence alleles in locus,
    1: one alternate allele on either strand in a locus, 
    2: two alternate allele on either strand in a locus (either as homozygous or compound heterozygous)
    '''
    
    # setup booleans
    mt = annotate_phased_entries(mt)
    homozygous = mt.a_homo
    compound_heterozygous = mt.a0_alt & mt.a1_alt
    heterozygous = mt.a0_alt | mt.a1_alt

    # aggregate onto gene level
    ht = (
        mt 
        .group_rows_by(mt.vep[gene_field])
        .aggregate(ko = hl.if_else( hl.agg.any(compound_heterozygous & ~homozygous) , 'CH', 
                        hl.if_else( hl.agg.any(homozygous & ~compound_heterozygous) , 'HO',
                        hl.if_else( hl.agg.any(compound_heterozygous & homozygous) , 'CH+HO',
                        hl.if_else( hl.agg.any(heterozygous) , 'HE', '')))))
        .filter_entries(mt.entries != '')
    )
    
    return ht

def extract_gene_ko_rows(burden_mt):
    assert ko in list(burden_mt)
    entries = burden_mt.entries()
    entries = entries.filter(~hl.literal('').contains(entries.ko))
    return entries

def extract_knockout_samples(mt, gene_field ='Gene', keep = ['HO','CH+HO']):
    r''' Collapses variants into genes for each individual, and returns
    a with three columns containing:
    1: The gene that was used for collapsing
    2: The current individual 
    4: Either CH (Compound heterozygous), HO (Homozygous) or CH+HO or H (Heterozygous).
       To avoid too long matrices, this is specified by the keep column.
    3: The variant and corresponding genotype for the variant in the sample,
      e.g. "chr22_50604850_G_A=1|1,chr22_50606762_C_T=1|1".
    '''
    
    # Parameters for KO/CH
    mt = annotate_phased_entries(mt)
    homozygous = mt.a_homo
    compound_heterozygous = mt.a0_alt & mt.a1_alt
    heterozygous = mt.a0_alt | mt.a1_alt

    # aggregate the gene-level
    ht1 = (mt 
          .group_rows_by(mt.vep[gene_field])
          .aggregate(ko = hl.if_else( hl.agg.any(compound_heterozygous & ~homozygous) , 'CH', 
                          hl.if_else( hl.agg.any(homozygous & ~compound_heterozygous) , 'HO',
                          hl.if_else( hl.agg.any(compound_heterozygous & homozygous) , 'CH+HO',
                          hl.if_else( hl.agg.any(heterozygous) , 'HE', '')))))
          .filter_entries(mt.entries != ''))

    # generate rsid to gt entries
    mt = mt.annotate_entries(rsid_entry = mt.rsid)
    mt = mt.annotate_entries(rsid_gt = hl.delimit([mt.rsid_entry, hl.str(mt.GT)], '='))

    # annotate only ko entries
    ht2 = (
            mt 
            .group_rows_by(mt.vep[gene_field])
            .aggregate(rsid = hl.agg.filter(mt.a_homo | (mt.a0_alt & mt.a1_alt), hl.agg.collect(mt.rsid_gt)))
    )

    # combine entries
    ent1 = ht1.entries()
    ent2 = ht2.entries()
    combined = ent1.annotate(genotypes = ent2[ent1[gene_field], ent1.s])
    combined = combined.filter(hl.set(keep).contains(combined.ko))
    line_merge = hl.delimit(combined.genotypes.rsid, ',')
    combined = combined.annotate(genotypes = combined.genotypes.annotate(rsid = line_merge))
    return combined


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
    '''
    Annotates entries with P(Knockout). Requires, that fields are
    already annotated with "singletons" count, and "dosage". 
    '''
    ko_mt = mt.annotate_entries(
        pKO = hl.if_else(
            mt.dosage == 2, 1, # knockout
            hl.if_else(
                mt.dosage == 1, 
                hl.if_else(mt.singletons >= 1, 1 - (1/2)**mt.singletons, 0), # one phased hetz
                hl.if_else(mt.singletons >= 2, 1 - 2*(1/2)**mt.singletons, 0), # zero phased hetz
            )
        )
    )
    return ko_mt

def get_prob_ko_matrix(mt_phased, mt_unphased, fields_drop = ['dosage','sigletons']):
    
    '''
    Annotates entries with P(Knockout). Requires a phased matrix and an unphased
    matrix that only contains singletons.
    '''
    
    # setup variables
    mt1 = mt_phased # contains phased non-singletons
    mt2 = mt_unphased # contains unphased singletons
    
    # add check to ensure that mt1 is phased
    # add check to ensure that mt2 only contains singletons
    
    # Determine probability of being KO given singletons and phased hetz
    mt1_burden = construct_phased_dosage_mt(mt1)
    mt2_burden = gene_burden_annotations_per_sample(mt2)
    mt_ko = mt1_burden.annotate_entries(singletons = mt2_burden[(mt1_burden.Gene, mt1_burden.s)].n)
    mt_ko = mt_ko.annotate_entries(singletons = hl.if_else(~hl.is_missing(mt_ko.singletons), mt_ko.singletons, 0 ))
    mt_ko = calc_p_ko(mt_ko)

    # drop not needed rows
    mt_ko = mt_ko.drop(*[f for f in fields_drop if f in mt_ko.entry])
    #mt_ko_entries = mt_ko.entries()
    #mt_ko_entries = mt_ko_entries.filter(~hl.is_missing(mt_ko_entries.pKO))
    return mt_ko

def get_dummy_by_dp(mt1, mt2, chrom):

    # mt1 is phased
    # mt2 is unphased

    # get probability matrix
    pmt = get_prob_ko_matrix(mt1, mt2, ["dosage","singletons"])
    pmt = pmt.annotate_entries(DP = pmt.pKO*2) # multiply probability by 2 (dosage encoded [0:2])
    pmt = pmt.drop('pKO')

    # create fake loci
    pmt = pmt.annotate_rows(locus = hl.parse_locus('chr' + str(chrom) + ':1'))
    pmt = pmt.annotate_rows(alleles = hl.literal(['X','Y']))
    pmt = pmt.annotate_rows(rsid = pmt.Gene)
    pmt = pmt.key_rows_by(pmt.locus, pmt.alleles)
    pmt = pmt.drop('Gene')
    return pmt
