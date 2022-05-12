import hail as hl
import argparse
import random

from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import variants
from ukb_utils import samples
from ko_utils import io

def main(args):

    chroms = args.chroms
    use_markers_by_kinship = args.use_markers_by_kinship
    out_prefix = args.out_prefix
    out_type = args.out_type
    use_markers_by_mac = args.use_markers_by_mac
    final_sample_list = args.final_sample_list
    sex = args.sex

    hail_init.hail_bmrc_init_local('logs/hail/_create_grm.log', 'GRCh37')
    hl._set_flags(no_whole_stage_codegen='1')
    mt = genotypes.get_ukb_genotypes_bed(chroms)
    mt = samples.remove_withdrawn(mt)
    markers = list() # keep track of final markers
    random.seed(1336) # seed for random

    if sex:
        mt = samples.filter_to_sex(mt, sex)

    if final_sample_list:
         ht_final_samples = hl.import_table(final_sample_list, no_header=True, key='f0',delimiter = ',')
         mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))

    if use_markers_by_mac:
        ht = genotypes.get_ukb_parsed_imputed_v3_mfi(chroms)
        mt_mac = mt.annotate_rows(info = ht[mt.row_key].info)
        mt_mac = mt_mac.annotate_rows(MAC = variants.get_mac_expr(mt_mac))
        mt_mac = mt_mac.filter_rows( (mt_mac.MAC >= 5)  & (mt_mac.MAC <= 20) & (mt_mac.info >= 0.5))
        mac_markers = mt_mac.rsid.collect()
        min_markers = min(int(use_markers_by_mac), len(mac_markers)) 
        print(f"(new) min_markers={min_markers}")
        markers.extend(random.sample(mac_markers, min_markers))

    if use_markers_by_kinship:
        ht = hl.import_table('/well/lindgren/UKBIOBANK/DATA/QC/ukb_snp_qc.txt', impute = True, delimiter = ' ')
        ht = ht.filter(ht.in_Relatedness == 1)
        kinship_markers = ht.rs_id.collect()
        markers.extend(kinship_markers) 

    mt = mt.filter_rows(hl.literal(set(markers)).contains(mt.rsid))
    io.export_table(mt, out_prefix, out_type)

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--chroms', nargs='+', default=None, help='chroms to be included')
    parser.add_argument('--use_markers_by_kinship', action='store_true', help='Use markers that have been used in UKBB as kinship markers')
    parser.add_argument('--use_markers_by_mac', default=None, help='Use N markers that have a MAC <= 20')
    parser.add_argument('--add_rare_variants', action='store_true', help='Add rare variants to plink file (this is required for SAIGE-GENE+ analysis')
    parser.add_argument('--final_sample_list', default=None, help='Path to HailTable that contains the final samples included in the analysis.')
    parser.add_argument('--sex', default=None, help='Only include "females" or "males".')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset (plink format)')
    parser.add_argument('--out_type', default=None, help='Path prefix for output dataset (plink format)')
    args = parser.parse_args()

    main(args)

