#!/usr/bin/env python3

#' @description Compound heterozygos HAIL pipeline
#' @todo Integrate SAIGE-GENE+
#' @todo Convert MatrixTable (sample, variant)-pairs to 'long' format
#' @todo Long format table should also contain variants found and their strand.
#' @DONE @todo Test burden_dosage function
#' @DONE @todo Create function that subsets variants of MODERATE impact
#' @todo hardy eq test for each variant
#' @todo add mt.repartition to different commands
#' @todo check worst consequence

import hail as hl
import os
import argparse

#WD = '/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb'
#DATA_DIR = f'{WD}/data'

def hail_init(chrom=None, log_prefix='get_vcf'):
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
        mt = hl.import_vcf(input_path, force_bgz=True, array_elements_required=False, min_partitions=50)
    elif input_type=='plink':
        mt = hl.import_plink(*[f'{input_path}.{x}' for x in ['bed','bim','fam']])
    if input_type!='mt' and cache:
        mt = mt.cache()
    return mt

def filter_max_maf(mt, maf=None):
    r'''Filter to variants to have maf less than {maf}'''
    mt = mt.filter_rows(hl.min(mt.info.AF)<maf)
    return mt

def filter_min_maf(mt, maf=None):
    r'''Filter to variants to have maf gt {maf}'''
    mt = mt.filter_rows(hl.min(mt.info.AF)>maf)
    return mt


def annotate_vep(mt, vep_path):
    r'''Annotate matrix table with VEP consequence from external file.'''
    print(f'Annotating with VEP file: {vep_path}')
    with open('data/vep/vep_fields.txt', 'r') as file:
        fields = file.read().strip().split(',')
    ht = hl.import_vcf(vep_path).rename({'info':'vep'}) 
    #ht = hl.annotate_rows(vep = hl.annoate(CSQ=ht.info.CSQ[0]))
    for i in range(len(fields)):
        ht = ht.annotate_rows(
            vep=ht.vep.annotate(
                col=ht.vep.CSQ.map(lambda x: (x.split('\\|')[i]))[0]
                ).rename({'col':f'{fields[i]}'})
        )
    # Extract various categories used downstream
    ht = ht.annotate_rows(vep = ht.vep.annotate(sift_pred = ht.vep.SIFT_pred.split('&')[0]))
    ht = ht.annotate_rows(vep = ht.vep.annotate(polyphen2_hdiv_pred = ht.vep.Polyphen2_HDIV_pred.split('&')[0]))
    ht = ht.annotate_rows(vep = ht.vep.annotate(polyphen2_hvar_pred = ht.vep.Polyphen2_HVAR_pred.split('&')[0]))

    #ht = ht.annotate_rows(vep = hl.struct(pred = ht.vep.SIFT_pred.split('&')[0]))
    #ht = ht.annotate(polyphen2 = hl.vep.struct(hdiv_pred = ht.vep.Polyphen2_HDIV_pred.split('&')[0]))
    #ht = ht.vep.annotate(polyphen2 = hl.vep.struct(hvar_pred = ht.vep.Polyphen2_HVAR_pred..split('&')[0]))
    # Categories for variants
    ptv = hl.set(["transcript_ablation", "splice_acceptor_variant",
              "splice_donor_variant", "stop_gained", "frameshift_variant"])
    missense = hl.set(["stop_lost", "start_lost", "transcript_amplification",
                   "inframe_insertion", "inframe_deletion", "missense_variant",
                   "protein_altering_variant", "splice_region_variant"])
    synonymous = hl.set(["incomplete_terminal_codon_variant", "stop_retained_variant", "synonymous_variant"])
    non_coding = hl.set(["coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant",
              "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant"])
    # create cases
    ht = ht.annotate_rows(consequence_category = 
        hl.case().when(ptv.contains(ht.vep.Consequence), "ptv")
             .when(missense.contains(ht.vep.Consequence) & 
                   (ht.vep.polyphen2_hdiv_pred == "") & 
                   (ht.vep.sift_pred == ""), "other_missense")
             .when(missense.contains(ht.vep.Consequence) & 
                   (ht.vep.polyphen2_hdiv_pred == "D") & 
                   (ht.vep.sift_pred == "D"), "damaging_missense")
             .when(missense.contains(ht.vep.Consequence), "other_missense")
             .when(synonymous.contains(ht.vep.Consequence), "synonymous")
             .when(non_coding.contains(ht.vep.Consequence), "non_coding")
             .default("NA")
    )

    #ht.filter_rows(hl.literal('chr22_16964821_G_A;chr22_16964821_G_C').contains(ht.rsid))

    
    # combine with matrix table
    mt = mt.annotate_rows(
        vep = ht.index_rows(mt.locus, mt.alleles).vep.drop('CSQ')
        )
    return(mt)



def filter_vep(mt, field, conditions):
    r'''Filter VEP field by condition(s) '''
    assert isinstance(conditions, list):
    assert field in list(mt.vep):
    conds = [[cond] for cond in conditions]
    mt = mt.filter_rows(hl.literal(conds).contains(mt.vep[field]))
    return mt

def write_sites(mt, out_prefix, keep_fields = []):
    r'''Write sites-only tsv
    `keep_fields` is a list of fields to keep in addition to chrom, pos, ref, alt
    '''
    ht = mt.rows().key_by()
    ht = ht.select(chrom = ht.locus.contig,
                   pos = ht.locus.position,
                   ref = ht.alleles[0],
                   alt = ht.alleles[1],
                   *keep_fields)
    ht.export(f'{out_prefix}.tsv')

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
    burden_mt = (
        mt 
        .group_rows_by(mt.vep[gene_field])
        .aggregate(dosage = hl.if_else( hl.agg.any((mt.a0_alt & mt.a1_alt) | mt.a_homo) , 2, 
                            hl.if_else( hl.agg.any((mt.a0_alt | mt.a1_alt)), 1, 0 )))
        #.burden_mt.filter_entries(burden_mt.entries != 0)
    )
    return burden_mt

def construct_summary_mt(mt, gene_field = 'Gene'):
    r''' Returns matrix table that contains dosage information from phased geneotypes.
    0: two refererence alleles in locus,
    1: one alternate allele on either strand in a locus, 
    2: two alternate allele on either strand in a locus (either as homozygous or compound heterozygous)
    '''
    mt = annotate_phased_entries(mt)
    ht = (
        mt 
        .group_rows_by(mt.vep[gene_field])
        .aggregate(dosage = hl.if_else( hl.agg.any((mt.a0_alt & mt.a1_alt) & (mt.a_homo == hl.literal(False))) , 4,  # compound hetz but not homo
                            hl.if_else( hl.agg.any((mt.a0_alt & mt.a1_alt) & (mt.a_homo == hl.literal(True))) , 3, # homozygous but not compoundd hetz
                            hl.if_else( hl.agg.any((mt.a0_alt | mt.a1_alt)) , 1, # heterozygous for either a1 or a2
                            hl.if_else( hl.agg.any((mt.a_homo == hl.literal(True))) , 2, 0))))) # homozygous
        .filter_entries(burden_mt.entries != 0)
    )    
    aggr = ht.aggregate_entries(hl.agg.counter(ht.dosage))
    print(aggr)
    return ht

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
    r'''Filter to samples in duos/trios, as defined by fam file
    '''
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


def main(args):
    
    # parser
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    out_type   = args.out_type
    vep_path   = args.vep_path
    
    chrom      = args.chrom
    max_maf    = args.max_maf
    min_maf    = args.min_maf 
    get_related = args.get_related
    get_unrelated = args.get_unrelated
    map_samples = args.map_samples
    get_europeans = args.get_europeans
    
    vep_impact = args.vep_impact
    vep_variant= args.vep_variant
    vep_loftee = args.vep_loftee
    vep_sites_write = args.vep_sites_write
    
    # run parser
    hail_init(chrom)
    mt = get_table(input_path=input_path, input_type=input_type)

    if max_maf:
        mt = filter_max_maf(mt, max_maf)

    if min_maf:
        mt = filter_min_maf(mt, min_maf)

    if get_related and not get_unrelated:
        mt = filter_to_unrelated(mt, get_related = True)

    if get_unrelated and not get_related:
        mt = filter_to_unrelated(mt, get_related = False)

    if vep_path:
        mt = annotate_vep(mt, vep_path)

    if vep_impact:
        mt = filter_vep('impact', vep_impact)

    if vep_variant:
        mt = filter_variants('variant', vep_variant)

    if vep_loftee:
        mt = filter_variants('loftee', vep_loftee)

    if map_samples:
        mt = translate_sample_ids(mt, 12788, 11867)

    if get_europeans:
        if map_samples is True:
            mt = get_genetically_european(mt)
        else:
            raise ValueError('EIDs are invalid. Did you map them using --map_samples?')

    if vep_sites_write:
        write_sites(mt=mt,
                    out_prefix=out_prefix,
                    keep_fields='info')

    if out_prefix:
        write(mt=mt,
              out_prefix=out_prefix,
              out_type=out_type)

    # test pipeline
    chrom=22
    mt = get_table('data/phased/ukb_wes_200k_phased_chr22.1of1.vcf.gz','vcf')
    #mt = filter_max_maf(mt, 0.02)
    mt = annotate_vep(mt, 'data/vep/output/ukb_wes_200k_vep_chr22.vcf')
    mt = filter_vep(mt, 'IMPACT', ['HIGH'])
    ht = construct_summary_mt(mt)

    #mt_burden = construct_phased_dosage_mt(mt)
    
    # sample filtering
    mt = filter_to_unrelated(mt, get_related = False)
    mt = translate_sample_ids(mt, 12788, 11867)
    mt = filter_to_european(mt)

    # generate LONG format described here: https://discuss.hail.is/t/how-to-write-a-matrixtable-to-a-file-as-a-tab-separated-table-in-wide-format/1338

    # How many individuals are compound hetz?
    
    mt = hl.import_vcf('data/tmp/ukb_wes_200k_phased_tmp_chr22.1of1.vcf.gz', force_bgz=True, array_elements_required=False, min_partitions=1)
    test = hl.vep(mt, 'data/vep/vep_env.json')

    mt = hl.import_vcf('data/phased/ukb_wes_200k_phased_chr22.1of1.vcf.gz', force_bgz=True, array_elements_required=False, min_partitions=50)
    test = hl.vep(mt, 'data/vep/vep_newest.json')
    mt = mt.filter_rows(hl.literal('chr22_16964821_G_A;chr22_16964821_G_C').contains(mt.rsid))
    mt = mt.repartition(5)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # I/O
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    # filtering variants
    parser.add_argument('--chrom', default=None, help='Chromosome to be used')
    parser.add_argument('--maf_max', default=None, help='Select all variants with a maf less than the indicated value')
    parser.add_argument('--min_maf', default=None, help='Select all variants with a maf greater than the indicated values')
    # filtering samples
    parser.add_argument('--get_related', default=None, help='Select all samples that are related')
    parser.add_argument('--get_unrelated', default=None, help='Select all samples that are unrelated')
    parser.add_argument('--map_samples', default=None, help='Map samples to lindgren UKBB app ID (11867)?')
    parser.add_argument('--get_europeans', default=None, help='Filter to genetically confimed europeans?')
    # VEP
    parser.add_argument('--vep_path', default=None, help='path to a .vcf file containing annotated entries by locus and alleles')
    parser.add_argument('--vep_impact', default=None, help='subset by VEP impact')
    parser.add_argument('--vep_variant', default=None, help='subset by VEP variant type (e.g. "stop_gained")')
    parser.add_argument('--vep_loftee', default=None, help='subset by VEP loftee flag, only HC or LC.')
    parser.add_argument('--vep_sites_write', default=None, help='Writes locus and allele information alongisde VEP annotations')



    args = parser.parse_args()

    main(args)


