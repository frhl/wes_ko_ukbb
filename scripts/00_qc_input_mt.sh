#!/usr/bin/env bash
#
#$ -N final_qc_input_mt
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/final_qc_input_mt.log
#$ -e logs/final_qc_input_mt.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 20
#$ -q short.qc
#$ -t 1


source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir_phased="/well/lindgren/UKBIOBANK/nbaya/wes_200k/phase_ukb_wes/data/phased/non_singleton"
readonly in_dir_unphased="data/unphased/post-qc"
readonly vep_dir="data/vep/full"
readonly gnomad_dir="/well/lindgren/flassen/ressources/gnomad/gnomad_v2_liftover/exomes"
readonly imputed_dir="/well/lindgren/UKBIOBANK/flassen/projects/ukb_compare/data/imputed/GRCh38"
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/qc_final"

# hail script
readonly hail_script="scripts/00_qc_input_mt.py"

# input paths
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_phased="${in_dir_phased}/ukb_wes_phased_non_singleton_chr${chr}-24xlong.qc-v4.2.2.vcf.gz"
readonly in_unphased="${in_dir_unphased}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly vep="${vep_dir}/ukb_wes_200k_full_vep_chr${chr}.vcf"
readonly gnomad="${gnomad_dir}/gnomad.exomes.r2.1.1.sites.${chr}.liftover_grch38.vcf.bgz"
readonly imputed="${imputed_dir}/ukb_imp_liftover_chr${chr}.mt"

# get final samples / variants from Duncan
readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'
readonly final_variant_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list'

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_chr${chr}"

# run hail
set_up_hail
set_up_vep
set_up_pythonpath_legacy  
python3 "${hail_script}" \
     --chrom ${chr} \
     --input_phased_path ${in_phased}\
     --input_unphased_path ${in_unphased} \
     --input_phased_type "vcf" \
     --input_unphased_type "mt" \
     --input_gnomad_path ${gnomad}\
     --input_imputed_path ${imputed}\
     --vep_path ${vep} \
     --final_sample_list ${final_sample_list} \
     --final_variant_list ${final_variant_list}\
     --out_prefix ${out_prefix}

print_update "Finished running HAIL for chr${chr}" "${SECONDS}"





