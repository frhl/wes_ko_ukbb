#!/usr/bin/env python3

import hail as hl
import os



def get_table(input_path, input_type, calc_info=True, cache=False):
    r'''Import mt/vcf/plink tables '''
    if input_type == 'mt':
        mt = hl.read_matrix_table(input_path)
    elif input_type == 'vcf':
        mt = hl.import_vcf(
            input_path,
            force_bgz=True,
            array_elements_required=False)
    elif input_type == 'plink':
        mt = hl.import_plink(
            *[f'{input_path}.{x}' for x in ['bed', 'bim', 'fam']])
    if input_type != 'mt' and cache:
        mt = mt.cache()
    if calc_info:
        mt = recalc_info(mt)
    return mt


def export_table(mt, out_prefix, out_type, checkpoint=False,
                 format_fields_to_drop=[]):
    r'''Merge file with VEP info '''
    assert out_type in {'mt', 'vcf', 'plink'}
    mt = mt.drop(*[f for f in format_fields_to_drop if f in mt.entry])
    if checkpoint:
        mt = mt.checkpoint(out_prefix + '.mt')
    if out_type == 'mt':
        print(f'\nWriting to MT: {out_prefix}.mt\n')
        mt.write(out_prefix + '.mt')
    elif out_type == 'vcf':
        out_vcf_path = out_prefix + '.vcf.bgz'
        print(f'\nExporting to VCF: {out_vcf_path}\n')
        if not hl.hadoop_is_file(out_vcf_path):
            metadata = get_vcf_metadata(
                format_fields_to_drop=format_fields_to_drop)
            hl.export_vcf(dataset=mt,
                          output=out_vcf_path,
                          parallel=None,
                          metadata=metadata,
                          tabix=False)
        else:
            raise ValueError(
                f'Cannot export VCF because {out_vcf_path} already exists')
    elif out_type == 'plink':
        print(f'\nExporting to PLINK: {out_prefix}.{{bed,bim,fam}}\n')
        mt = annotate_with_fam_fields(mt)
        if not any(hl.hadoop_is_file(out_prefix + s)
                   for s in ['bed', 'bim', 'fam']):  # don't overwrite if any files exist
            hl.export_plink(dataset=mt,
                            output=out_prefix,
                            fam_id=mt.fam_id,
                            pat_id=mt.pat_id,
                            mat_id=mt.mat_id,
                            is_female=mt.is_female,
                            varid=mt.rsid)
        else:
            raise ValueError(
                f'Cannot export to PLINK because at least one of {out_prefix}.{{bed,bim,fam}} already exists')


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


def translate_sample_ids(ht, from_app: int, to_app: int):
    r'''Translate sample IDs from one UKB application to another
    `from_app` and `to_app` are UKB application IDs. This function only supports
    translation in either direction between applications 11867 (Lindgren) and 12788 (McVean)
    '''
    from_app, to_app = int(from_app), int(to_app)
    valid_apps = {11867, 12788}
    assert (
        from_app in valid_apps) & (
        to_app in valid_apps), f'from_app and to_app must be in {valid_apps}'
    assert from_app != to_app, 'from_app and to_app must be different'
    print(f'Mapping IDs from UKB app {from_app} to {to_app}')
    id_dict = hl.import_table(
        '/well/lindgren/UKBIOBANK/nbaya/resources/ukb11867_to_ukb12788.sample_ids.txt',
        delimiter='\\s+',
        key=f'eid_{from_app}')
    ht = ht.key_cols_by(s=id_dict[ht.s][f'eid_{to_app}'])
    undefined_ct = ht.aggregate_cols(hl.agg.sum(hl.is_missing(ht.s)))
    assert undefined_ct == 0, f'[translate_sample_ids]: Not all sample IDs mapped perfectly ({undefined_ct}/{ht.count()} IDs are undefined)'
    return ht


def annotate_european(mt, genetically_european=True):
    r'''annotate white british (app 11867) /well/lindgren/UKBIOBANK/DATA/QC/ukb_sqc_v2.txt
    or genetically european by projecting ancestries into 1KG prpject data'''
    if genetically_european:
        ht1 = hl.import_table(
            '/well/lindgren/dpalmer/ukb_get_EUR/data/final_EUR_list.tsv',
            no_header=True).rename(
            {
                'f0': 'eid'}).key_by('eid')
        ht1 = ht1.annotate(eur=1)
        mt = mt.annotate_cols(eur=ht1[mt.s].eur)
    else:
        ht2 = hl.import_table('/well/lindgren/flassen/ressources/ukb/white_british/210921_ukbb_white_british_samples.txt',
                              types={'eid': hl.tstr, 'in.white.British.ancestry.subset': hl.tint32}).key_by('eid')
        mt = mt.annotate_cols(
            eur=ht2[mt.s]['in.white.British.ancestry.subset'])
    return mt


def filter_to_european(mt, genetically_european=True, use_existing_field=True):
    r'''Get white british (app 11867) /well/lindgren/UKBIOBANK/DATA/QC/ukb_sqc_v2.txt
    or genetically european by projecting ancestries into 1KG prpject data'''
    # if 'eur' not in list(mt.cols) or not use_existing_field:
    mt = annotate_european(mt, genetically_european)
    undefined_eur = mt.aggregate_cols(hl.agg.sum(hl.is_missing(mt.eur)))
    pre_filter_count = mt.count()
    if undefined_eur == pre_filter_count[1]:
        raise ValueError(
            '[get_european]: IDs for europeans does not match keys in MatrixTable!')
    if undefined_eur > 0:
        print(
            f'[get_european]: Not all samples IDs mapped perfectly ({undefined_eur}/{pre_filter_count[1]} IDs are undefined)')
    mt = mt.filter_cols(mt.eur == 1)
    post_filter_count = mt.count()
    print(
        f'[get_european]:{post_filter_count[1]}/{pre_filter_count[1]} IDs were included as genetically european.')

    return mt


def filter_to_european_legacy(
        mt, genetically_european=True, only_annotate=False):
    r'''Get white british (app 11867) /well/lindgren/UKBIOBANK/DATA/QC/ukb_sqc_v2.txt
    and genetically european from /well/lindgren/UKBIOBANK/laura/k_means_clustering_pcs/ukbb_genetically_european_k4_4PCs_self_rep_Nov2020.txt'''

    # filter to either genetically european or UKBB
    if genetically_european:
        ht1 = hl.import_table('/well/lindgren/UKBIOBANK/laura/k_means_clustering_pcs/ukbb_genetically_european_k4_4PCs_self_rep_Nov2020.txt',
                              types={'eid': hl.tstr, 'genetically_european': hl.tint32}).key_by('eid')
        mt = mt.annotate_cols(eur=ht1[mt.s].genetically_european)
    else:
        ht2 = hl.import_table('/well/lindgren/flassen/ressources/ukb/white_british/210921_ukbb_white_british_samples.txt',
                              types={'eid': hl.tstr, 'in.white.British.ancestry.subset': hl.tint32}).key_by('eid')
        mt = mt.annotate_cols(
            eur=ht2[mt.s]['in.white.British.ancestry.subset'])
    # count and subset
    undefined_eur = mt.aggregate_cols(hl.agg.sum(hl.is_missing(mt.eur)))
    pre_filter_count = mt.count()
    if undefined_eur == pre_filter_count[1]:
        raise ValueError(
            '[get_european]: IDs for europeans does not match keys in MatrixTable!')
    if undefined_eur > 0:
        print(
            f'[get_european]: Not all samples IDs mapped perfectly ({undefined_eur}/{pre_filter_count[1]} IDs are undefined)')
    if only_annotate == False:
        mt = mt.filter_cols(mt.eur == 1)
        post_filter_count = mt.count()
        print(
            f'[get_european]:{post_filter_count[1]}/{pre_filter_count[1]} IDs were included as genetically european.')
    return mt


def get_fam(app_id=12788, wes_200k_only=False):
    # use files created by the command:
    #   /well/lindgren/UKBIOBANK/nbaya/wes_200k/phase_ukb_wes/utils/make_fam.py --make_fam
    # github:
    # https://github.com/nikbaya/phase_ukb_wes/blob/a5a97de8604d4ea1c7e7cfeecd9a1f57c2af2357/utils/make_fam.py
    assert app_id in {11867, 12788}
    fam_path = f'/well/lindgren/UKBIOBANK/nbaya/resources/ukb{app_id}_{"wes_200k_" if wes_200k_only else ""}pedigree.fam'
    fam = hl.import_table(
        paths=fam_path,
        key='IID',
        types={
            f: 'str' for f in [
                'FID',
                'IID',
                'PAT',
                'MAT',
                'SEX',
                'PHEN']})
    return fam


def recalc_info(mt, maf=None, info_field='info', gt_field='GT'):
    r'''Recalculate INFO fields AF, AC, AN and keep sites with minor allele
    frequency > `maf The fields AF and AC are made into integers instead of
    an array, only showing the AF and AC for the alt allele
    '''
    if info_field in mt.row:
        mt = mt.annotate_rows(
            info=mt[info_field].annotate(**hl.agg.call_stats(mt[gt_field],
                                                             mt.alleles)
                                         ).drop('homozygote_count')
        )  # recalculate AF, AC, AN, drop homozygote_count
    else:
        mt = mt.annotate_rows(
            **{info_field: hl.agg.call_stats(mt[gt_field], mt.alleles)})
    mt = mt.annotate_rows(
        **{info_field: mt[info_field].annotate(
            **{field: mt[info_field][field][1] for field in ['AF', 'AC']}
        )}
    )  # keep only the 2nd entries in AF and AC, which correspond to the alt alleles
    if maf is not None:
        mt = mt.filter_rows(
            (mt[info_field].AF > maf) & (
                mt[info_field].AF < (
                    1 - maf)))
    return mt


def filter_to_unrelated(mt, get_related=False, maf=None):
    r'''Filter to samples in duos/trios, as defined by fam file'''
    fam = get_fam()
    # need to filter to samples present in mt first before collecting FIDs
    fam = fam.filter(hl.is_defined(mt.cols()[fam.IID]))
    fam_ct = fam.count()
    col_ct = mt.count_cols()
    # samples with missing IIDs cannot be included in the related sample set
    print(f'\n{col_ct-fam_ct}/{col_ct} samples have IIDs missing from fam file')
    iids = fam.IID.collect()
    offspring_fids = fam.filter(hl.literal(iids).contains(fam.PAT)
                                | hl.literal(iids).contains(fam.MAT)).FID.collect()  # get the FIDs of offspring with at least one parent in the dataset
    # subset to samples which share FIDs with offspring
    fam = fam.filter(hl.literal(offspring_fids).contains(fam.FID))
    if get_related:
        mt = mt.filter_cols(hl.is_defined(fam[mt.s]))  # related
    else:
        mt = mt.filter_cols(~hl.is_defined(fam[mt.s]))  # unrelated
    mt = recalc_info(mt=mt, maf=maf)
    return mt


def filter_hwe(mt, cut_off=1e-12):
    mt = mt.annotate_rows(hwe=hl.agg.hardy_weinberg_test(mt.GT))
    rows = mt.hwe.p_values >= cut_off
    mt = mt.filter_rows(rows)
    return mt


def get_vcf_metadata(info_fields_to_drop=[],
                     format_fields_to_drop=['RNC', 'PL']):
    r'''Merge file with VEP info '''
    path = '/well/ukbb-wes/pvcf/oqfe/ukbb-wes-oqfe-pvcf-chr1.vcf.gz'  # use a VCF that has the same fields and is unlikely to be deleted
    metadata = hl.get_vcf_metadata(path)
    for f in info_fields_to_drop:
        metadata['info'].pop(f)
    for f in format_fields_to_drop:
        metadata['format'].pop(f)
    return metadata


def annotate_with_fam_fields(mt):
    r'''Annotate `mt` with FID (fam_id), PAT (pat_id), MAT (mat_id), SEX (is_female)from fam file with all UKB samples'''
    fam = get_fam(wes_200k_only=False)
    mt = mt.annotate_cols(fam_id=fam[mt.s].FID,
                          pat_id=fam[mt.s].PAT,
                          mat_id=fam[mt.s].MAT,
                          is_female=hl.if_else(fam[mt.s].SEX == '2',
                                               True,
                                               hl.if_else(fam[mt.s].SEX == '1',
                                                          False,
                                                          hl.null(hl.tbool))))
    return mt


def annotate_snpid(mt, delim='_'):
    r'''Annotate field snpid based on current locus and alleles'''
    ids = hl.delimit([mt.locus.contig, hl.str(
        mt.locus.position), mt.alleles[0], mt.alleles[1]], delim)
    mt = mt.annotate_rows(snpid=ids)
    return mt


def annotate_rsid(
        mt, dbsnp_path='/well/lindgren/flassen/ressources/dbsnp/GRCh38/155/GCF_000001405.39.gz', build='GRCh38'):
    r'''Use dbSNP to annotate all rsIDs in the a matrix table.'''
    recode = {f"NC_0000{i}.{j}": f"chr{i}" for i in (
        list(range(1, 23)) + ['X', 'Y']) for j in ('09', '10', '11', '12', '13', '14')}
    dbsnp = hl.import_vcf(dbsnp_path,
                          reference_genome=build,
                          contig_recoding=recode,
                          skip_invalid_loci=True,
                          force_bgz=True)

    rsids = dbsnp.index_rows(mt.locus, mt.alleles).rsid
    mt = mt.annotate_rows(rsid=rsids)
    return mt


def default_to_snpid_when_missing_rsid(mt):
    r'''rsid is converted to snpid when it is missing'''
    return mt.annotate_rows(rsid=hl.if_else(
        hl.is_missing(mt.rsid), mt.snpid, mt.rsid))


def is_phased(mt):
    ''' Check if the input contains phased data. Returns Bool'''
    mt = mt.annotate_entries(phased=(mt.GT == hl.parse_call('0|0')) |
                             (mt.GT == hl.parse_call('1|0')) |
                             (mt.GT == hl.parse_call('0|1')) |
                             (mt.GT == hl.parse_call('1|1'))
                             )
    return mt.aggregate_entries(hl.agg.any(mt.phased))

def union_phased_with_unphased(mt_phased, mt_unphased):
    r'''Union rows of phased haplotypes with unphased singetons. 
    
    Union the rows of phased haplotypes with unphased singletons. Only
    overlapping entries will be kept. Only samples defined in mt_phased
    will be kept.
    
    :param mt_phased: a phased matrix table with info.AC row
    :param mt_unphased: a phased matrix_table with info row.
    '''
    
    # Filter to singletons
    mt_unphased = mt_unphased.filter_rows(mt_unphased.info.AC == 1)
    
    # drop rows not required
    mt_phased = mt_phased.drop('info')
    mt_unphased = mt_unphased.drop('info')
    
    # subset to overlapping entries
    overlapping_entries = list(set(mt_phased.entry) & set(mt_unphased.entry))
    mt_unphased = mt_unphased.select_entries(*overlapping_entries)
    mt_phased = mt_phased.select_entries(*overlapping_entries)
    
    # subset to overlapping samples (assuming unphased is larger)
    defined_in_phased = hl.is_defined(mt_phased.cols()[mt_unphased.s])
    mt_unphased = mt_unphased.filter_cols(defined_in_phased)
    
    # keep track of which variants were phased
    mt_unphased = mt_unphased.annotate_rows(phased=0)
    mt_phased = mt_phased.annotate_rows(phased=1)
    
    return mt_phased.union_rows(mt_unphased)


