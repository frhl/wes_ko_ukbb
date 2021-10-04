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
        hl.case().when(hl.literal(set(PLOF_CSQS)).contains(mt.vep.worst_csq_by_gene_canonical.most_severe_consequence), "ptv")
             .when(hl.literal(set(MISSENSE_CSQS)).contains(mt.vep.worst_csq_by_gene_canonical.most_severe_consequence) & 
                   (~hl.is_defined(mt.dbnsfp.cadd_phred_score) | 
                    ~hl.is_defined(mt.dbnsfp.revel_score)), "other_missense")                                   
             .when(hl.literal(set(MISSENSE_CSQS)).contains(mt.vep.worst_csq_by_gene_canonical.most_severe_consequence) & 
                   (mt.dbnsfp.cadd_phred_score >= 20) & 
                   (mt.dbnsfp.revel_score >= 0.6), "damaging_missense") 
             .when(hl.literal(set(MISSENSE_CSQS)).contains(mt.vep.worst_csq_by_gene_canonical.most_severe_consequence), "other_missense")
             .when(hl.literal(set(SYNONYMOUS_CSQS)).contains(mt.vep.worst_csq_by_gene_canonical.most_severe_consequence), "synonymous")
             .when(hl.literal(set(OTHER_CSQS)).contains(mt.vep.worst_csq_by_gene_canonical.most_severe_consequence), "non_coding")
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
    #in_mt = in_mt.explode_rows(in_mt.vep.worst_csq_by_gene_canonical)
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
    #in_mt = in_mt.explode_rows(in_mt.vep.worst_csq_by_gene_canonical)
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
    #mt = mt.explode_rows(mt.vep.worst_csq_by_gene_canonical)
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
    #mt = mt.explode_rows(mt.vep.worst_csq_by_gene_canonical)
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
    #mt = mt.explode_rows(mt.vep.worst_csq_by_gene_canonical)
    mt = mt.group_rows_by(
        gene_id = mt.vep.worst_csq_by_gene_canonical.gene_id
        ).aggregate(n = hl.agg.count_where(mt.GT.is_non_ref()))
    return mt

def gene_burden_category_annotations_per_sample(mt):
    r''' calculate gene burden by counting variants in gene stratified by variant category'''
    #mt = mt.explode_rows(mt.vep.worst_csq_by_gene_canonical)
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
    #mt = mt.explode_rows(mt.vep.worst_csq_by_gene_canonical)
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

(base) [mmq446@rescomp1 ko_utils]$ ls
__init__.py  __init__.pyc  __pycache__  analysis.py  legacy.py  phenos.py  phenos.pyc  qc.py
(base) [mmq446@rescomp1 ko_utils]$ cat qc.py 
#!/usr/bin/env python3

import hail as hl
import os

def hail_init(chrom=None, log_prefix='wes_analysis'):
    r'''Initialize Hail '''
    assert chrom in range(1,23), 'only autosomes allowed'
    n_slots = os.environ.get('NSLOTS', 1)
    chr_suffix = '' if chr is None else f'_chr{chrom}'
    WD = '/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb'
    hl.init(log=f'{WD}/logs/hail-{log_prefix}{chr_suffix}.log',
            default_reference='GRCh38',
            master=f'local[{n_slots}]')

def get_table(input_path, input_type, cache=False):
    r'''Import mt/vcf/plink tables '''
    if input_type=='mt':
        mt = hl.read_matrix_table(input_path)
    elif input_type=='vcf':
        mt = hl.import_vcf(input_path, force_bgz=True, array_elements_required=False, min_partitions=100)
    elif input_type=='plink':
        mt = hl.import_plink(*[f'{input_path}.{x}' for x in ['bed','bim','fam']])
    if input_type!='mt' and cache:
        mt = mt.cache()
    mt = recalc_info(mt)
    return mt

def export_table(mt, out_prefix, out_type, checkpoint=False, format_fields_to_drop=[]):
    r'''Merge file with VEP info '''
    assert out_type in {'mt', 'vcf', 'plink'}
    mt = mt.drop(*[f for f in format_fields_to_drop if f in mt.entry])
    if checkpoint:
        mt = mt.checkpoint(out_prefix+'.mt')
    if out_type=='mt':
        print(f'\nWriting to MT: {out_prefix}.mt\n')
        mt.write(out_prefix+'.mt')
    elif out_type=='vcf':
        out_vcf_path = out_prefix+'.vcf.bgz'
        print(f'\nExporting to VCF: {out_vcf_path}\n')
        if not hl.hadoop_is_file(out_vcf_path):
            metadata = get_vcf_metadata(format_fields_to_drop=format_fields_to_drop)
            hl.export_vcf(dataset=mt,
                          output=out_vcf_path,
                          parallel=None,
                          metadata = metadata,
                          tabix = False)
        else:
            raise ValueError(f'Cannot export VCF because {out_vcf_path} already exists')
    elif out_type=='plink':
        print(f'\nExporting to PLINK: {out_prefix}.{{bed,bim,fam}}\n')
        mt = annotate_with_fam_fields(mt)
        if not any(hl.hadoop_is_file(out_prefix+s) for s in ['bed','bim','fam']): # don't overwrite if any files exist
            hl.export_plink(dataset=mt,
                            output=out_prefix,
                            fam_id=mt.fam_id,
                            pat_id=mt.pat_id,
                            mat_id=mt.mat_id,
                            is_female=mt.is_female,
                            varid=mt.rsid)
        else:
            raise ValueError(f'Cannot export to PLINK because at least one of {out_prefix}.{{bed,bim,fam}} already exists')



def filter_max_af(mt, maf=None):
    r'''Filter to variants to have af less than {maf}'''
    mt = mt.filter_rows(hl.min(mt.info.AF)<af)
    return mt

def filter_min_af(mt, maf=None):
    r'''Filter to variants to have af gt {maf}'''
    mt = mt.filter_rows(hl.min(mt.info.AF)>af)
    return mt

def filter_max_maf(mt, maf=None):
    r'''boolean for variants that have maf less than {maf}'''
    maf_expr = (mt.info.AF<maf) | (mt.info.AF>(1-maf))
    return maf_expr

def filter_min_maf(mt, maf=None):
    r'''boolean for variants that have maf gt {maf}'''
    maf_expr = (mt.info.AF>maf) & (mt.info.AF<(1-maf))
    return maf_expr

def filter_maf(mt, min_maf = None, max_maf = None):
    r'''Filter to variants based on a certain min/max MAF threshold'''
    if min_maf is not None:
            mt = mt.filter_rows(filter_min_maf(mt, min_maf))
    if max_maf is not None:
            mt = mt.filter_rows(filter_max_maf(mt, max_maf))
    return mt

def filter_max_mac(mt, mac=None):
    r'''Filter to variants to have maf less than {maf}'''
    if mac is not None:
        mt = mt.filter_rows(hl.min(mt.info.AC)<=mac)
    return mt

def filter_min_mac(mt, mac=None):
    r'''Filter to variants to have mac gt {mac}'''
    if mac is not None:
        mt = mt.filter_rows(hl.min(mt.info.AC)>mac)
    return mt

def filter_min_missing(mt, value = 0.05):
    r'''Filter variants to have less than {value} in genotype missigness"'''
    assert value >= 0 and value < 1
    missing = hl.agg.mean(hl.is_missing(mt.GT)) <= value
    mt = mt.filter_rows(missing)
    return mt

def translate_sample_ids(ht, from_app: int, to_app: int):
    r'''Translate sample IDs from one UKB application to another
    `from_app` and `to_app` are UKB application IDs. This function only supports
    translation in either direction between applications 11867 (Lindgren) and 12788 (McVean)
    '''
    from_app, to_app = int(from_app), int(to_app)
    valid_apps = {11867, 12788}
    assert (from_app in valid_apps) & (to_app in valid_apps), f'from_app and to_app must be in {valid_apps}'
    assert from_app!=to_app, 'from_app and to_app must be different'
    print(f'Mapping IDs from UKB app {from_app} to {to_app}')
    id_dict = hl.import_table('/well/lindgren/UKBIOBANK/nbaya/resources/ukb11867_to_ukb12788.sample_ids.txt',delimiter='\s+', key=f'eid_{from_app}')
    ht = ht.key_cols_by(s = id_dict[ht.s][f'eid_{to_app}'])
    undefined_ct = ht.aggregate_cols(hl.agg.sum(hl.is_missing(ht.s)))
    assert undefined_ct==0, f'[translate_sample_ids]: Not all sample IDs mapped perfectly ({undefined_ct}/{ht.count()} IDs are undefined)'
    return ht

def filter_to_european(mt, genetically_european = True, only_annotate = False):
    r'''Get white british (app 11867) /well/lindgren/UKBIOBANK/DATA/QC/ukb_sqc_v2.txt
    and genetically european from /well/lindgren/UKBIOBANK/laura/k_means_clustering_pcs/ukbb_genetically_european_k4_4PCs_self_rep_Nov2020.txt''' 
    
    # filter to either genetically european or UKBB 
    if genetically_european:
        ht1 = hl.import_table('/well/lindgren/UKBIOBANK/laura/k_means_clustering_pcs/ukbb_genetically_european_k4_4PCs_self_rep_Nov2020.txt',
            types={'eid': hl.tstr, 'genetically_european': hl.tint32}).key_by('eid')
        mt = mt.annotate_cols(eur = ht1[mt.s].genetically_european)
    else:
        ht2 = hl.import_table('/well/lindgren/flassen/ressources/ukb/white_british/210921_ukbb_white_british_samples.txt',
            types={'eid': hl.tstr, 'in.white.British.ancestry.subset': hl.tint32}).key_by('eid')
        mt = mt.annotate_cols(eur = ht2[mt.s]['in.white.British.ancestry.subset'])
    # count and subset
    undefined_eur = mt.aggregate_cols(hl.agg.sum(hl.is_missing(mt.eur)))
    pre_filter_count = mt.count()
    if undefined_eur == pre_filter_count[1]:
        raise ValueError('[get_european]: IDs for europeans does not match keys in MatrixTable!')
    if undefined_eur > 0:
        print(f'[get_european]: Not all samples IDs mapped perfectly ({undefined_eur}/{pre_filter_count[1]} IDs are undefined)')
    if only_annotate == False:
        mt = mt.filter_cols(mt.eur == 1)
        post_filter_count = mt.count()
        print(f'[get_european]:{post_filter_count[1]}/{pre_filter_count[1]} IDs were included as genetically european.')
    return mt

def get_fam(app_id=12788, wes_200k_only=False):
    # use files created by the command:
    #   /well/lindgren/UKBIOBANK/nbaya/wes_200k/phase_ukb_wes/utils/make_fam.py --make_fam
    # github: https://github.com/nikbaya/phase_ukb_wes/blob/a5a97de8604d4ea1c7e7cfeecd9a1f57c2af2357/utils/make_fam.py
    assert app_id in {11867,12788}
    fam_path = f'/well/lindgren/UKBIOBANK/nbaya/resources/ukb{app_id}_{"wes_200k_" if wes_200k_only else ""}pedigree.fam'
    fam = hl.import_table(paths=fam_path, key='IID',types={f:'str' for f in ['FID','IID','PAT','MAT','SEX','PHEN']})
    return fam

def recalc_info(mt, maf=None, info_field='info', gt_field='GT'):
    r'''Recalculate INFO fields AF, AC, AN and keep sites with minor allele
    frequency > `maf The fields AF and AC are made into integers instead of 
    an array, only showing the AF and AC for the alt allele
    '''
    if info_field in mt.row:
        mt = mt.annotate_rows(
            info = mt[info_field].annotate(**hl.agg.call_stats(mt[gt_field],
                                                               mt.alleles)
                                               ).drop('homozygote_count')
            ) # recalculate AF, AC, AN, drop homozygote_count
    else:
        mt = mt.annotate_rows(**{info_field: hl.agg.call_stats(mt[gt_field], mt.alleles)})
    mt = mt.annotate_rows(
        **{info_field: mt[info_field].annotate(
            **{field:mt[info_field][field][1] for field in ['AF','AC']}
            )}
        ) # keep only the 2nd entries in AF and AC, which correspond to the alt alleles
    if maf is not None:
        mt = mt.filter_rows((mt[info_field].AF>maf) & (mt[info_field].AF<(1-maf)))
    return mt

def filter_to_unrelated(mt, get_related=False, maf=None):
    r'''Filter to samples in duos/trios, as defined by fam file'''
    fam = get_fam()
    fam = fam.filter(hl.is_defined(mt.cols()[fam.IID])) # need to filter to samples present in mt first before collecting FIDs
    fam_ct = fam.count()
    col_ct = mt.count_cols()
    print(f'\n{col_ct-fam_ct}/{col_ct} samples have IIDs missing from fam file') # samples with missing IIDs cannot be included in the related sample set
    iids = fam.IID.collect()
    offspring_fids = fam.filter(hl.literal(iids).contains(fam.PAT)
                                |hl.literal(iids).contains(fam.MAT)).FID.collect() # get the FIDs of offspring with at least one parent in the dataset
    fam = fam.filter(hl.literal(offspring_fids).contains(fam.FID)) # subset to samples which share FIDs with offspring
    if get_related:
        mt = mt.filter_cols(hl.is_defined(fam[mt.s])) # related
    else:
        mt = mt.filter_cols(~hl.is_defined(fam[mt.s])) # unrelated
    mt = recalc_info(mt=mt, maf=maf)
    return mt

def filter_hwe(mt, cut_off = 1e-12):
    mt =  mt.annotate_rows(hwe = hl.agg.hardy_weinberg_test(mt.GT))
    rows = mt.hwe.p_values >= cut_off
    mt = mt.filter_rows(rows)
    return mt

def get_vcf_metadata(info_fields_to_drop=[], format_fields_to_drop=['RNC','PL']):
    r'''Merge file with VEP info '''
    path = '/well/ukbb-wes/pvcf/oqfe/ukbb-wes-oqfe-pvcf-chr1.vcf.gz' # use a VCF that has the same fields and is unlikely to be deleted
    metadata = hl.get_vcf_metadata(path)
    for f in info_fields_to_drop:
        metadata['info'].pop(f)
    for f in format_fields_to_drop:
        metadata['format'].pop(f)
    return metadata

def annotate_with_fam_fields(mt):
    r'''Annotate `mt` with FID (fam_id), PAT (pat_id), MAT (mat_id), SEX (is_female)from fam file with all UKB samples'''
    fam = get_fam(wes_200k_only=False)
    mt = mt.annotate_cols(fam_id = fam[mt.s].FID,
                          pat_id = fam[mt.s].PAT,
                          mat_id = fam[mt.s].MAT,
                          is_female = hl.if_else(fam[mt.s].SEX=='2',
                                                 True,
                                                 hl.if_else(fam[mt.s].SEX=='1',
                                                            False,
                                                            hl.null(hl.tbool))))
    return mt


def annotate_snpid(mt, delim = '_'):
    r'''Annotate field snpid based on current locus and alleles'''
    ids = hl.delimit([mt.locus.contig, hl.str(mt.locus.position), mt.alleles[0], mt.alleles[1]], delim)
    mt = mt.annotate_rows(snpid = ids)
    return mt

def annotate_rsid(mt, dbsnp_path = '/well/lindgren/flassen/ressources/dbsnp/GRCh38/155/GCF_000001405.39.gz', build = 'GRCh38'):
    r'''Use dbSNP to annotate all rsIDs in the a matrix table.'''
    recode = {f"NC_0000{i}.{j}":f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y']) for j in ('09','10','11','12','13','14')}
    dbsnp = hl.import_vcf(dbsnp_path, 
                       reference_genome=build, 
                       contig_recoding=recode, 
                       skip_invalid_loci=True,
                       force_bgz=True)
    
    rsids = dbsnp.index_rows(mt.locus, mt.alleles).rsid
    mt = mt.annotate_rows(rsid = rsids)
    return mt


def default_to_snpid_when_missing_rsid(mt):
    r'''rsid is converted to snpid when it is missing'''
    return mt.annotate_rows(rsid = hl.if_else(hl.is_missing(mt.rsid), mt.snpid, mt.rsid))

def is_phased(mt):
    ''' Check if the input contains phased data. Returns Bool'''
    mt = mt.annotate_entries(phased = 
                        (mt.GT ==  hl.parse_call('0|0')) |
                        (mt.GT ==  hl.parse_call('1|0')) |
                        (mt.GT ==  hl.parse_call('0|1')) |
                        (mt.GT ==  hl.parse_call('1|1'))
                       )
    return mt.aggregate_entries(hl.agg.any(mt.phased))



