import os
import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import samples
from ko_utils import qc

def main(args):

    chroms = args.chroms
    subset_markers_by_kinship = args.subset_markers_by_kinship
    out_prefix = args.out_prefix
    final_sample_list = args.final_sample_list
    sex = args.sex

    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh37')
    hl._set_flags(no_whole_stage_codegen='1')
    mt = genotypes.get_ukb_genotypes_bed(chroms)
    mt = samples.remove_withdrawn(mt)

    if sex:
        mt = samples.filter_to_sex(mt, sex)

    if final_sample_list:
         ht_final_samples = hl.import_table(final_sample_list, no_header=True, key='f0',delimiter = ',')
         mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))
    
    hl.export_plink(mt, out_prefix, ind_id = mt.s)
   
if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chroms', nargs='+', default=None, help='chroms to be included')
    parser.add_argument('--subset_markers_by_kinship', action='store_true', help='Only use markers that have been used in UKBB as kinship markers')
    parser.add_argument('--final_sample_list', default=None, help='Path to HailTable that contains the final samples included in the analysis.')
    parser.add_argument('--sex', default=None, help='Only include "females" or "males".')
    parser.add_argument('--phenotypes', default=None, help='Path to phenotype table that also includes sex.')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset (plink format)')
    args = parser.parse_args()

    main(args)

