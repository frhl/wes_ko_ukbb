#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import samples
from ko_utils import qc
from ko_utils import analysis

def main(args):
    
    # parser
    input_annotation_path = args.input_annotation_path
    out_prefix = args.out_prefix
    final_variant_list = args.final_variant_list

    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    
    ht = hl.read_table(input_annotation_path)
    
    if final_variant_list:
        ht_final_variants = hl.import_table(final_variant_list, types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
        ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)
        ht = ht.filter(hl.is_defined(ht_final_variants[ht.key]))
    
    gene_map_path = out_prefix + '_gene_map.ht'
    gene_map_ht = analysis.create_gene_map_ht(ht)
    gene_map_ht.write(gene_map_path, overwrite=True)

    gene_map_processed_path = out_prefix + "_gene_map_processed.ht"
    gene_map_ht = hl.read_table(gene_map_path)
    gene_map_ht = analysis.post_process_gene_map_ht(gene_map_ht)
    gene_map_ht.write(gene_map_processed_path, overwrite=True)   
    gene_map_ht = hl.read_table(gene_map_processed_path)

    groups="ptv,damaging_missense|ptv_LC,ptv|damaging_missense|ptv_LC,ptv|damaging_missense,damaging_missense,other_missense,synonymous"
    # Create a distinct file for each annotation
    for group in groups.split(','):
        lst = group.split('|')
        #gene_ht = gene_map_ht.filter(group==gene_map_ht.annotation)
        gene_ht = gene_map_ht.fiter(hl.literal(set(lst)).contains(gene_map_ht.annotation))
        output_group_input_path_tmp = out_prefix + '_' + group.replace('|','_') + '.tsv.gz'
        gene_ht.select(
            group=gene_ht.gene_id + '_' + gene_ht.gene_symbol + '_' + gene_ht.annotation,
            variant=hl.delimit(gene_ht.variants, '\t')
        ).key_by().drop('start').export(output_group_input_path_tmp, header=False)

    
    output_vep_ht_path = out_prefix + '_vep.ht'
    output_summary_path = out_prefix + '_vep_summary.ht'
    analysis.count_variants(output_vep_ht_path, vep_vcf_path)

    ht = hl.read_table(output_vep_ht_path)
    ht = process_consequences(ht)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    ht = ht.group_by(
        gene=ht.vep.worst_csq_by_gene_canonical.gene_symbol,
        consequence=annotation_case_builder(ht.vep.worst_csq_by_gene_canonical)
        ).partition_hint(100).aggregate(n_variants=hl.agg.count())
    ht.write(output_summary_path, overwrite=True)




if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--input_annotation_path', default=None, help='path to HailTable with VEP and dbNSFP annotations')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    parser.add_argument('--final_variant_list', default=None, help='Path to hail table with final variants to be included')     

    args = parser.parse_args()

    main(args)

