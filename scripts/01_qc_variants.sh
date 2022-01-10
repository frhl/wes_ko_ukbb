#!/usr/bin/env bash
#
#$ -N qc_union_variants
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/qc_union_variants.log
#$ -e logs/qc_union_variants.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q long.qc@@long.hge
#$ -t 1-22

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir_phased="/well/lindgren/UKBIOBANK/nbaya/wes_200k/phase_ukb_wes/data/phased/non_singleton"
readonly in_dir_unphased="data/unphased/post-qc"
readonly gnomad_dir="/well/lindgren/flassen/ressources/gnomad/gnomad_v2_liftover/exomes"
readonly imputed_dir="/well/lindgren/UKBIOBANK/flassen/projects/ukb_compare/data/imputed/GRCh38"
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/variants"
readonly chr=$( get_chr ${SGE_TASK_ID} ) 

# hail script
readonly hail_script="scripts/01_qc_variants.py"

# input tables
readonly in_phased="${in_dir_phased}/ukb_wes_phased_non_singleton_chr${chr}-24xlong.qc-v4.2.2.vcf.gz"
readonly in_unphased="${in_dir_unphased}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly gnomad="${gnomad_dir}/gnomad.exomes.r2.1.1.sites.${chr}.liftover_grch38.vcf.bgz"
readonly imputed="${imputed_dir}/ukb_imp_liftover_chr${chr}.mt"
readonly annotation_table="data/vep/hail/ukb_wes_200k_chr${chr}_vep.ht"
readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'
readonly final_variant_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list'

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_chr${chr}"

# run hail
mkdir -p ${out_dir}
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
     --final_sample_list ${final_sample_list} \
     --final_variant_list ${final_variant_list}\
     --out_prefix ${out_prefix}\
     --combine_datasets

print_update "Finished running HAIL for chr${chr}" "${SECONDS}"





