#!/usr/bin/env bash
#
#$ -N create_mt
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/create_mt.log
#$ -e logs/create_mt.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 3
#$ -q long.qc@@long.hga
#$ -t 1-22

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly in_dir_phased="/well/lindgren/UKBIOBANK/nbaya/wes_200k/phase_ukb_wes/data/phased/non_singleton"
readonly in_dir_unphased="data/unphased/post-qc"
readonly spark_dir="data/tmp/spark_dir"
readonly out_dir="data/mt"

readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_phased="${in_dir_phased}/ukb_wes_phased_non_singleton_chr${chr}-24xlong.qc-v4.2.2.vcf.gz"
readonly in_unphased="${in_dir_unphased}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly annotation_table="data/vep/hail/ukb_wes_200k_chr${chr}_vep.ht"

readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'
readonly final_variant_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list'

readonly out_prefix="${out_dir}/ukb_wes_200k_merged_chr${chr}"
readonly out="${out_prefix}.mt"

readonly hail_script="scripts/03_create_mt.py"

mkdir -p ${out_dir}
if [ ! -f "${out}/_SUCCESS" ]; then
  set_up_hail
  set_up_pythonpath_legacy  
  warn_spark_dir_size
  python3 "${hail_script}" \
     --chrom ${chr} \
     --input_phased_path ${in_phased}\
     --input_unphased_path ${in_unphased} \
     --input_phased_type "vcf" \
     --input_unphased_type "mt" \
     --input_annotation_path ${annotation_table}\
     --final_sample_list ${final_sample_list} \
     --final_variant_list ${final_variant_list}\
     --out_prefix ${out_prefix} \
     --dbsnp "155" \
     --varid \
     --out_type "mt"
else
  print_update "file ${out} already exists. Skipping!"
fi



print_update "Finished running HAIL for chr${chr}" "${SECONDS}"





