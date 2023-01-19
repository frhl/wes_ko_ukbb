#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import io
from ko_utils import ko


def main(args):

    input_path = args.input_path
    input_type = args.input_type
    clinvar_path = args.clinvar_path
    out_prefix = args.out_prefix

    # import phased/unphased data
    hail_init.hail_bmrc_init(log='logs/hail/knockout.log', default_reference='GRCh38', min_block_size=128)
       
    # load clinvar
    clinvar_path = "/well/lindgren/flassen/ressources/clinvar/ftp//clinvar_20230107.vcf.gz"
    recoding = {str(i):"chr" + str(i) for i in range(1,23)}
    mt1 = hl.import_vcf(clinvar_path, contig_recoding=recoding, force=True, skip_invalid_loci=True)
    ht1 = mt1.rows()

    # load data to annotate 
    mt2 = io.import_table(input_path, input_type, calc_info=False) 
    ht2 = mt2.rows()
    ht2 = ht2.select(*[ht2.info, ht2.consequence.vep.worst_csq_by_gene_canonical, ht2.consequence_category])
    ht2 = ht2.annotate(clinvar=ht1[ht2.key].info)
        
    ht2.flatten().export(out_prefix + ".txt.gz")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--input_path', default=None, help='Path to input')
    parser.add_argument('--input_type', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--clinvar_path', default=None, help='Input type, either "mt", "vcf" or "plink"')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')

    args = parser.parse_args()

    main(args)



