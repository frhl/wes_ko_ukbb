#!/usr/bin/env python3

#' @description Compound heterozygos HAIL pipeline
#' @todo Integrate SAIGE-GENE+
#' @DONE @todo Convert MatrixTable (sample, variant)-pairs to 'long' format
#' @DONE @todo Long format table should also contain variants found and their strand.
#' @DONE @todo Test burden_dosage function
#' @DONE @todo Create function that subsets variants of MODERATE impact
#' @DONE @todo hardy eq test for each variant
#' @DONE @todo add mt.repartition to different commands
#' @DONE @todo check worst consequence
#' @DONE @todo export ultra rare variants
#' @todo check which allele is included in VEP? What about MAF thresholds?


import hail as hl
import argparse
import pandas
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
                          metadata = metadata)
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

## QC PIPELINE
def filter_max_maf(mt, maf=None):
    r'''Filter to variants to have maf less than {maf}'''
    mt = mt.filter_rows(hl.min(mt.info.AF)<maf)
    return mt

def filter_min_maf(mt, maf=None):
    r'''Filter to variants to have maf gt {maf}'''
    mt = mt.filter_rows(hl.min(mt.info.AF)>maf)
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

def filter_to_european(mt):
    r'''Get white british (app 11867) /well/lindgren/UKBIOBANK/DATA/QC/ukb_sqc_v2.txt
    and genetically european from /well/lindgren/UKBIOBANK/laura/k_means_clustering_pcs/ukbb_genetically_european_k4_4PCs_self_rep_Nov2020.txt''' 
    ht = hl.import_table('/well/lindgren/UKBIOBANK/laura/k_means_clustering_pcs/ukbb_genetically_european_k4_4PCs_self_rep_Nov2020.txt',
        types={'eid': hl.tstr, 'genetically_european': hl.tint32}).key_by('eid')
    mt = mt.annotate_cols(eur = ht[mt.s].genetically_european)
    undefined_eur = mt.aggregate_cols(hl.agg.sum(hl.is_missing(mt.eur)))
    pre_filter_count = mt.count()
    if undefined_eur == pre_filter_count[1]:
        raise ValueError('[get_genetically_european]: IDs for genetically europeans does not match keys in MatrixTable!')
    if undefined_eur > 0:
        print(f'[get_genetically_european]: Not all samples IDs mapped perfectly ({undefined_eur}/{pre_filter_count[1]} IDs are undefined)')
    #mt = mt.filter_cols(mt.s[mt.eur])
    mt = mt.filter_cols(mt.eur == 1)
    post_filter_count = mt.count()
    print(f'[get_genetically_european]:{post_filter_count[1]}/{pre_filter_count[1]} IDs were included as genetically european.')
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



def remake_rsid(mt, delim = '_'):
    r'''Re-create rsids based on current locus and alleles'''
    ids = hl.delimit([mt.locus.contig, hl.str(mt.locus.position), mt.alleles[0], mt.alleles[1]], delim)
    mt = mt.annotate_rows(rsid = ids)
    return mt

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

def calc_ko_prob(mt):
    '''
    Calculate the probability of a gene being a KO in an invidual given the 
    amount of phased hetz in a gene (construct_phased_dosage_mt) and the singletons 
    in the gene (gene_burden_annotations_per_sample). 
    '''
    #z = x.annotate_entries(k = y[(x.Gene, x.s)].n)
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



def main(args):
    
    # parser
    input_phased_path = args.input_phased_path
    input_phased_type = args.input_phased_type
    input_unphased_path = args.input_unphased_path
    input_unphased_type = args.input_unphased_type
    out_prefix = args.out_prefix
    out_type   = args.out_type
    vep_path   = args.vep_path
    
    chrom      = int(args.chrom)
    maf_max    = (args.maf_max)
    maf_min    = (args.maf_min)
    hwe        = (args.hwe)
    missing    = args.missing
    
    get_related = args.get_related
    get_unrelated = args.get_unrelated
    get_europeans = args.get_europeans
    vep_variants = args.vep_variants
    ko_matrix = args.ko_matrix
    ko_samples = args.ko_samples
    export_burden = args.export_burden

    # run parser
    hail_init(chrom)
    mt1 = get_table(input_path=input_phased_path, input_type=input_phased_type) # 12788
    mt2 = get_table(input_path=input_unphased_path, input_type=input_unphased_type) # 11867 (for singletons)

    ### Sample filtering

    if get_related and not get_unrelated:
        mt1 = filter_to_unrelated(mt1, get_related = True)
        mt2 = filter_to_unrelated(mt2, get_related = True)

    if get_unrelated and not get_related:
        mt1 = filter_to_unrelated(mt1, get_related = False)
        mt2 = filter_to_unrelated(mt2, get_related = False)

    if get_europeans:
        mt1 = translate_sample_ids(mt1, 12788, 11867)
        mt1 = filter_to_european(mt1)
        mt2 = filter_to_european(mt2)

    ### Variant filtering/annotations
    # Using mt2 as a singleton refereence, so remove those with AC > 1
    mt2 = filter_max_mac(mt2, 1)

    if missing:
        mt1 = filter_min_missing(mt1, 0.05)
        mt2 = filter_min_missing(mt2, 0.05)

    if maf_max:
        mt1 = filter_max_maf(mt1, float(maf_max))

    if maf_min:
        mt1 = filter_min_maf(mt1, float(maf_min))

    if hwe:
        mt1 = filter_hwe(mt1, float(hwe))

    if vep_path:
        mt1 = annotate_vep(mt1, vep_path)
        mt2 = annotate_vep(mt2, vep_path)

    #### get stats

    if export_burden:

        # Count burden per gene per individual
        mt1_cat = gene_burden_category_annotations_per_sample(mt1)
        mt2_cat = gene_burden_category_annotations_per_sample(mt2)

        # combine singleton table and full table
        res = mt1_cat.annotate_entries(singletons = mt2_cat[(mt1_cat.Gene, mt1_cat.consequence_category), mt1_cat.s].n)
        res = res.annotate_entries(singletons = hl.if_else(hl.is_missing(res.singletons),0,res.singletons))
        res = res.annotate_entries(total = res.n + hl.if_else(hl.is_missing(res.singletons),0,res.singletons))
        res = res.entries()
        res = res.filter_rows(res.total > 0)

        # export data
        res.export(out_prefix + '_burden.tsv.bgz')

    
    #if vep_variants:
    #    summary_ptv = summarize_variants(mt, 'ptv')
    #    summary_missense = summarize_variants(mt, 'damaging_missense')
    #    df = pd.DataFrame([summary_ptv, summary_missense], \
    #        columns=['Singletons','non-singletons','total AC', 'Genes'],\
    #        index = ['ptv','missense'])
    #    df.to_csv(out_prefix + '_plof_variants.csv', index=True)

    #if ko_samples:
    #    mt_ko_sample = extract_knockout_samples(mt)
    #    mt_ko_sample.export(prefix + '_ko_samples.tsv.bgz')

    #if ko_matrix:
    #    mt_ko_matrix = construct_phased_dosage_mt(mt)
    #    mt_ko_matrix.export(prefix + '_ko_matrx.tsv.bgz')

    #if out_prefix & out_type:
    #    export_table(mt, out_prefix, out_type)



if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_phased_path', default=None, help='Path to input')
    parser.add_argument('--input_phased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--input_unphased_path', default=None, help='Path to input that contains singletons')
    parser.add_argument('--input_unphased_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    # filtering variants
    parser.add_argument('--chrom', default=None, help='Chromosome to be used')
    parser.add_argument('--maf_max', default=None, help='Select all variants with a maf less than the indicated value')
    parser.add_argument('--maf_min', default=None, help='Select all variants with a maf greater than the indicated values')
    parser.add_argument('--hwe', default=None, help='Filter variants by HWE threshold')
    parser.add_argument('--missing', default=0.05, help='Filter variants by missingness threshold')
    # filtering samples
    parser.add_argument('--get_related', action='store_true', help='Select all samples that are related')
    parser.add_argument('--get_unrelated', action='store_true', help='Select all samples that are unrelated')
    parser.add_argument('--get_europeans', action='store_true', help='Filter to genetically confimed europeans?')
    # out
    parser.add_argument('--export_burden', action='store_true', help='Export burden variant count by gene and and individuals.')
    parser.add_argument('--vep_path', default=None, help='path to a .vcf file containing annotated entries by locus and alleles')
    parser.add_argument('--vep_variants', action='store_true', help='Generate a summary of filter variants')
    parser.add_argument('--ko_samples', action='store_true', help='Get the genes/individuals that are KO and the SNPs involved')
    parser.add_argument('--ko_matrix', action='store_true', help='Generate a gene x sample matrix with KO status')
    
    args = parser.parse_args()

    main(args)



