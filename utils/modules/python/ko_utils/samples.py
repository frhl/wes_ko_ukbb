#!/usr/bin/env python3

import hail as hl


def get_phenotype_path():
    r''' returns the current path to Duncan's Phenotype table'''
    return("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/phenotypes/curated_covar_phenotypes_cts.tsv.gz")

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


def filter_to_females(mt):
    '''Filter MatrixTable with chrX variants
    :param mt: MatrixTable to filter
    '''
    # Sex is genetic sex (UKB field 22001), unless:
    #   - Genetic sex is missing, in which case reported sex (UKB field 31) is used.
    #   - If genetic and reported sex are defined and genetic sex does not
    #     match reported sex, then sex is set to unknown/missing.
    # Sex has been recoded to PLINK coding (2=female, 1=male, 0=unknown)
    # from UKB coding (0=female, 1=male).
    ht = hl.import_table(
            get_phenotype_path(), 
            impute = True, 
            key = 'eid', missing = ["NA",""], 
            types = {"eid": hl.tstr},
            force = True
            )
    # Filter to females
    return mt.filter_cols(ht[mt.s].sex == 0)


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




