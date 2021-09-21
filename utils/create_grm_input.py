import os
import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes
from ko_utils import qc

def main(args):

    hail_init.hail_bmrc_init('logs/hail_grm.log', 'GRCh37')
    
    # input variables
    chroms = args.chroms
    subset_markers_by_kinship = args.subset_markers_by_kinship
    out_prefix = args.out_prefix
    subset_samples_by_eur = args.subset_samples_by_eur
    subset_samples_by_wes200k = args.subset_samples_by_wes200k   
     
    # combine multiple matrix tables
    mt = genotypes.get_ukb_imputed_v3_bgen(chroms)
       
    # createse a sparse GRM
    if subset_markers_by_kinship:
        ht = hl.import_table('/well/lindgren/UKBIOBANK/DATA/QC/ukb_snp_qc.txt', impute = True, delimiter = ' ')
        ht = ht.filter(ht.in_Relatedness == 1)
        rsids = ht.rs_id.collect()
        mt = mt.filter_rows(hl.literal(rsids).contains(mt.rsid))

    # subset population by WES data samples
    if subset_samples_by_wes200k:
        ids = genotypes.get_ukb_wes_200k_post_qc_path(chrom=22).s.collect()
        mt = mt.filter_cols(hl.literal(ids).contains(mt.s))
    
    # get white british ancestry
    if subset_samples_by_eur:
        mt = qc.filter_to_european(mt)
    
    hl.export_plink(mt, out_prefix, ind_id = mt.s)
    
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chroms', nargs='+', default=None, help='chroms to be included')
    parser.add_argument('--subset_markers_by_kinship', action='store_true', help='Only use markers that have been used in UKBB as kinship markers')
    parser.add_argument('--subset_samples_by_eur', action='store_true', help='Subset samples to those that are genetically european.') 
    parser.add_argument('--subset_samples_by_wes200k', default=None, help='Subset samples to those that are presennt in the post QC WES.')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset (plink format)')
    args = parser.parse_args()

    main(args)

