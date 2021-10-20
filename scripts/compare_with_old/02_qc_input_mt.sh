#!/usr/bin/env bash
#
#$ -N old_qc_input_mt
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/old_qc_input_mt.log
#$ -e logs/old_qc_input_mt.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qc


source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir_phased="/gpfs3/well/lindgren/UKBIOBANK/nbaya/wes_200k/phase_ukb_wes/data/old-all/all_variants/phased/24slots-lindgren.qe/1000000-1000-1666666/chr21"
readonly in_dir_unphased="/gpfs3/well/lindgren/UKBIOBANK/nbaya/wes_200k/phase_ukb_wes/data/old-all/all_variants/unphased"
readonly gnomad_dir="/well/lindgren/flassen/ressources/gnomad/gnomad_v2_liftover/exomes"
readonly imputed_dir="/well/lindgren/UKBIOBANK/flassen/projects/ukb_compare/data/imputed/GRCh38"
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/debug_tmp"

# hail script
readonly hail_script="scripts/02_qc_input_mt.py"

# input paths
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_phased="${in_dir_phased}/ukb_wes_200k_phased_chr1.1of1.vcf.gz"
readonly in_unphased="${in_dir_unphased}/ukb_wes_hail_chr*.vcf.gz"
readonly gnomad="${gnomad_dir}/gnomad.exomes.r2.1.1.sites.21.liftover_grch38.vcf.bgz"
readonly imputed="${imputed_dir}/ukb_imp_liftover_chr21.mt"
readonly annotation_table="data/debug_tmp/ukb_wes_200k_chr21_vep.ht"

# get final samples / variants from Duncan
readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'
readonly final_variant_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list'

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_chr21"

# run hail
set_up_hail
set_up_pythonpath_legacy  
python3 "${hail_script}" \
     --chrom ${chr} \
     --input_phased_path ${in_phased}\
     --input_unphased_path ${in_unphased} \
     --input_phased_type "vcf" \
     --input_unphased_type "mt" \
     --input_gnomad_path ${gnomad}\
     --input_imputed_path ${imputed}\
     --input_annotation_path ${annotation_table} \
     --out_prefix ${out_prefix}
     #--final_sample_list ${final_sample_list} \
     #--final_variant_list ${final_variant_list}\

print_update "Finished running HAIL for chr${chr}" "${SECONDS}"





