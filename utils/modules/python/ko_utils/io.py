#!/usr/bin/env python3

import hail as hl


def import_table(input_path, input_type, calc_info=True, cache=False):
    r'''Import mt/vcf/plink tables '''
    assert input_type in ['mt','vcf','plink']
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

def rbind_matrix_tables(mt1: hl.MatrixTable, mt2: hl.MatrixTable, 
                verbose: bool = True):
    r'''rowbind two matrix tables. 

    Note: if the samples are not ordered you will need
    to call ukb_utils::tables::order_cols on the MatrixTables
    before running this function.
    
    :param mt1: A MatrixTable (e.g. with phased entries)
    :param mt2: A MatrixTable (e.g. with unphased entries)
    '''
    
    
    mt1 = mt1.drop('info')
    mt2 = mt2.drop('info')

    overlap = set(mt1.row) & set(mt2.row)
    overlap = overlap - set(list(mt1.row_key))
    mt1 = mt1.select_rows(*overlap)
    mt2 = mt2.select_rows(*overlap)
    if verbose:
        print(f"Note: {overlap} row field(s) kept.")
     
    overlap = list(set(mt1.entry) & set(mt2.entry))
    mt1 = mt1.select_entries(*overlap)
    mt2 = mt2.select_entries(*overlap)
    if verbose:
        print(f"Note: {overlap} entry fields kept.")

    s1 = mt1.s.collect()
    s2 = mt2.s.collect()
    overlap = list(set(s1) & set(s2))
    mt1 = mt1.filter_cols(hl.literal(overlap).contains(mt1.s))
    mt2 = mt2.filter_cols(hl.literal(overlap).contains(mt2.s))
    if verbose:
        n = len(overlap)
        n_max = max(len(s1), len(s2))
        pct = (n / n_max)*100
        print(f"Note: {n}/{n_max} ({pct:.2f}%) overlapping samples kept.")

    return mt1.union_rows(mt2)




