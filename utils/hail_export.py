#!/usr/bin/env python3

import hail as hl
import os
import argparse

WD = '/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb'
DATA_DIR = f'{WD}/data'

def hail_init(chrom=None, log_prefix='get_vcf'):
    r'''Initialize Hail '''
    n_slots = os.environ.get('NSLOTS', 1)
    chr_suffix = '' if chr is None else f'_chr{chrom}'
    WD = '/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb'
    hl.init(log=f'{WD}/logs/hail-{log_prefix}{chr_suffix}.log',
            default_reference='GRCh38',
            master=f'local[{n_slots}]')

def import_table(input_path, input_type, cache=False):
    r'''Import mt/vcf/plink tables '''
    if input_type=='mt':
        mt = hl.read_matrix_table(input_path)
    elif input_type=='vcf':
        mt = hl.import_vcf(input_path, force_bgz=True, array_elements_required=False)
    elif input_type=='plink':
        mt = hl.import_plink(*[f'{input_path}.{x}' for x in ['bed','bim','fam']])
    if input_type!='mt' and cache:
        mt = mt.cache()
    return mt

def filter_min_maf(mt, maf=None):
    r'''Filter to variants with minor allele frequency '''
    if maf is not None:
        if mt.info.AF.dtype==hl.dtype('array<float64>'):
            mt = mt.filter_rows(hl.max(mt.info.AF)<maf)
        else:
            raise ValueError('MatrixTable does not have an info.AF field with dtype array<float64>')
    return mt

def filter_min_maf(mt, maf=None):
    r'''Filter to variants with minor allele frequency '''
    if maf is not None:
        if mt.info.AF.dtype==hl.dtype('array<float64>'):
            mt = mt.filter_rows(hl.min(mt.info.AF)>maf)
        else:
            raise ValueError('MatrixTable does not have an info.AF field with dtype array<float64>')
    return mt


   #mt = mt.filter_rows(hl.agg.call_stats(mt.GT, mt.alleles).AF[1] > maf)
#def filter_samples(mt, subset = 'WB'):
#    r'''Filter samples to only include QCED individuals in reeference population '''
#    if subset=='WB':
#        qt = hl.import_table('/well/lindgren/UKBIOBANK/DATA/WES_200k/ukb23155_c1_b0_v1_s200632.fam')
#        # still need to select the WB samples?
#        qt = qt.rename({'eid' : 's'}).key_by('s')
#        mt = mt.filter_cols(hl.is_defined(qt[mt['s']]))
#        #post_filter_count = mt.count()
#        return(mt)
#    else:
#        raise TypeError('Subset is not valid. Must be "WB"')

def annotate_vep(mt):
    r'''Merge file with VEP info '''
    ht = hl.import_vcf('derived/vep/output/ukb_wes_200k_vep_chr22.vcf')
    ht = ht.annotate_rows(info=ht.info.annotate(ensgid=ht.info.CSQ.map(lambda x: x.split('\\|')[0])))
    ht = ht.annotate_rows(info=ht.info.annotate(enstid=ht.info.CSQ.map(lambda x: x.split('\\|')[1])))
    ht = ht.annotate_rows(info=ht.info.annotate(variant=ht.info.CSQ.map(lambda x: x.split('\\|')[3])))
    ht = ht.annotate_rows(info=ht.info.annotate(impact=ht.info.CSQ.map(lambda x: x.split('\\|')[4])))
    ht = ht.annotate_rows(info=ht.info.annotate(loftee_flag=ht.info.CSQ.map(lambda x: x.split('\\|')[15])))
    # merge with vep
    mt = mt.annotate_rows(info=mt.info.annotate(ensgid=ht.index_rows(mt.locus, mt.alleles).info.ensgid))
    mt = mt.annotate_rows(info=mt.info.annotate(variant=ht.index_rows(mt.locus, mt.alleles).info.variant))
    mt = mt.annotate_rows(info=mt.info.annotate(impact=ht.index_rows(mt.locus, mt.alleles).info.impact))
    mt = mt.annotate_rows(info=mt.info.annotate(loftee_flag=ht.index_rows(mt.locus, mt.alleles).info.loftee_flag))
    return(mt)

def select_field(mt, field = 'impact', condition = 'HIGH'):
    mt = mt.filter_rows(mt.info[field].contains(condition))
    return mt

def annotate_gt(mt):
    mt = mt.annotate_rows(info=mt.info.annotate(hl.agg.count_where(mt.GT.is_het())))
    mt = mt.annotate_rows(info=mt.info.annotate(hl.agg.count_where(mt.GT.is_homo())))
    # if vep.info.gene is present, annotate compound hetz
    return(mt)

def get_fam(app_id=11867, wes_200k_only=False):
    # use files created by the command:
    #   /well/lindgren/UKBIOBANK/nbaya/wes_200k/phase_ukb_wes/utils/make_fam.py --make_fam
    # github: https://github.com/nikbaya/phase_ukb_wes/blob/a5a97de8604d4ea1c7e7cfeecd9a1f57c2af2357/utils/make_fam.py
    assert app_id in {11867,12788}
    fam_path = f'/well/lindgren/UKBIOBANK/nbaya/resources/ukb{app_id}_{"wes_200k_" if wes_200k_only else ""}pedigree.fam'
    fam = hl.import_table(paths=fam_path, key='IID',types={f:'str' for f in ['FID','IID','PAT','MAT','SEX','PHEN']})
    return fam

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
    
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    out_type   = args.out_type
    chrom      = args.chrom
    max_maf    = args.max_maf
    min_maf    = args.min_maf 
    
    # filtering data
    chrom=22
    hail_init(chrom)
    mt = import_table('data/phased/ukb_wes_200k_phased_chr22.1of1.vcf.gz','vcf')
    mt = filter_maf(mt, 0.02, 'lower')
    mt = annotate_vep(mt)
    mt = select(mt, 'impact','HIGH')

    export_table
    # 
    # find compound hetz? phased data?
    # export plink file


    mt.show(width = 80)



    mt = read_input(input_path=input_path,
                    input_type=input_type)




if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # I/O
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None, help='Type of output dataset (options: mt, vcf, plink)')
    # filtering
    parser.add_argument('--chrom', default=None, help='Chromosome to be used')
    parser.add_argument('--maf_max', default=None, help='Select all variants with a maf less than the indicated value')
    parser.add_argument('--min_maf', default=None, help='Select all variants with a maf greater than the indicated values')

    args = parser.parse_args()

    main(args)