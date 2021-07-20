#!/usr/bin/env python3

import hail as hl
import os
import argparse
import re

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

def filter_maf(mt, maf=None, bounds = 'lower'):
    r'''Filter to variants with minor allele frequency '''
    if maf is not None:
        if mt.info.AF.dtype==hl.dtype('array<float64>'):
            if bounds=='lower':
                mt = mt.filter_rows(hl.min(mt.info.AF)<maf)
            elif bounds=='upper':
                mt = mt.filter_rows(hl.min(mt.info.AF)>maf)
            else:
                raise TypeError('Bounds must be either "upper" or "lower"')
        else:
            raise ValueError('MatrixTable does not have an info.AF field with dtype array<float64>')
    return mt


def filter_samples(mt, subset = 'WB'):
    r'''Filter samples to only include QCED individuals in reeference population '''
    if subset=='WB':
        qt = hl.import_table('/well/lindgren/UKBIOBANK/stefania/RelGroups/2021_03_12/QCWB.txt')
        qt = qt.rename({'eid' : 's'}).key_by('s')
        mt = mt.filter_cols(hl.is_defined(qt[mt['s']]))
        #post_filter_count = mt.count()
        return(mt)
    else:
        raise TypeError('Subset is not valid. Must be "WB"')

def join_with_vep(mt, vep_path, keyword = None):
    r'''Merge file with VEP info '''
    vep = hl.import_vcf('derived/vep/output/ukb_wes_200k_vep_chr22.vcf', array_elements_required=False)
    mt = mt.annotate_rows(info=mt.info.annotate(CSQ = vep.index_rows(mt.locus, mt.alleles).info.CSQ))
    # select rows by specific keyword, e.g. HIGH / MODERATE
    #mt = mt[]
    return(mt)

def rename_rsid():
    # rename rsids
    return None

def count_homozygous():
    return None

def export_table():
    # export bed files that need to be analyzed
    return None




    




def main(args):
    chrom      = args.chrom
    input_path = args.input_path
    input_type = args.input_type
    out_prefix = args.out_prefix
    out_type   = args.out_type

    # filtering data
    hail_init(chrom)
    mt = import_table('data/phased/ukb_wes_200k_phased_chr22.1of1.vcf.gz','vcf')
    mt = filter_maf(mt, 0.02, 'lower')
    mt = filter_samples(mt, 'WB')
    mt = join_with_vep(mt)



    mt = read_input(input_path=input_path,
                    input_type=input_type)




if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', default=None, help='Chromosome to filter')
    args = parser.parse_args()

    main(args)