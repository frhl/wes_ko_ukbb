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
    add_rare_variants = args.add_rare_variants
    final_sample_list = args.final_sample_list
    subset_samples_by_genet_eur = args.subset_samples_by_genet_eur
    subset_samples_by_ukbb_eur = args.subset_samples_by_ukbb_eur
    subset_samples_by_wes200k = args.subset_samples_by_wes200k   
     
    # combine multiple matrix tables
    mt = genotypes.get_ukb_imputed_v3_bgen(chroms)   
    
    # only keep final samples
    if final_sample_list:
         ht_final_samples = hl.import_table(final_sample_list, no_header=True, key='f0')
         mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))

    # createse a sparse GRM
    if subset_markers_by_kinship:
        ht = hl.import_table('/well/lindgren/UKBIOBANK/DATA/QC/ukb_snp_qc.txt', impute = True, delimiter = ' ')
        ht = ht.filter(ht.in_Relatedness == 1)
        rsids = ht.rs_id.collect()
        mt = mt.filter_rows(hl.literal(rsids).contains(mt.rsid))

    # subset population by WES data samples
    if subset_samples_by_wes200k:
        path = genotypes.get_ukb_wes_200k_post_qc_path(chrom=22)
        ids = qc.get_table(path, 'mt').s.collect()
        mt = mt.filter_cols(hl.literal(ids).contains(mt.s))
    
    # get white british ancestry (k-means)
    if subset_samples_by_genet_eur or subset_samples_by_ukbb_eur:
        mt = qc.filter_to_european(mt, genetically_european = subset_samples_by_genet_eur)
    
    n = mt.count()
    print(f"Count after subsetting: {n}")

    hl.export_plink(mt, out_prefix, ind_id = mt.s)
   
    # add rare variants for SAIGE 
    if add_rare_variants:
        
        # get high quality rare variants
        rsids = hl.literal('')
        for c in range(1,24):
            mfi = genotypes.get_ukb_imputed_v3_mfi(c)
            mfi = mfi.filter((mfi.maf < 0.001) & 
                     (mfi.maf > 0.000001) &
                     (mfi.info == 1))
            rsid = mfi.rsid
            rsids += rsid
        rsids = hl.set(rsids.collect())
        # get genotype file and filter
        mt2 = genotypes.get_ukb_imputed_v3_bgen(chroms, entry_fields = ['GT'])
        mt2 = mt2.filter_cols(hl.set(mt.s.collect()).contains(mt.s))
        mt2 = mt2.filter_rows(rsids.contains(mt2.rsid))
        n = mt2.count()
        print(f"added {n} rare variants.")
        mt_out = mt.union_rows(mt2)
        hl.export_plink(mt_out, out_prefix + "_rare", ind_id = mt_out.s) 


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chroms', nargs='+', default=None, help='chroms to be included')
    parser.add_argument('--subset_markers_by_kinship', action='store_true', help='Only use markers that have been used in UKBB as kinship markers')
    parser.add_argument('--subset_samples_by_wes200k', action='store_true', help='Subset samples to those that are presennt in the post QC WES.')
    parser.add_argument('--subset_samples_by_genet_eur', action='store_true', help='Subset samples to those that are genetically european.') 
    parser.add_argument('--subset_samples_by_ukbb_eur', action='store_true', default = False, help='Subset samples to those that are european designated by UKBB.') 
    parser.add_argument('--add_rare_variants', action='store_true', help='Add rare variants to plink file (this is required for SAIGE-GENE+ analysis')
    parser.add_argument('--final_sample_list', default=None, help='Path to HailTable that contains the final samples included in the analysis.')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset (plink format)')
    args = parser.parse_args()

    main(args)

