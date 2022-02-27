#!/usr/bin/env bash
#
#$ -N combine_mt
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/combine_mt.log
#$ -e logs/combine_mt.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q long.qc@@long.hga
#$ -t 21

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly in_dir_phased="/well/lindgren/UKBIOBANK/nbaya/wes_200k/phase_ukb_wes/data/phased/non_singleton"
readonly in_dir_unphased="data/unphased/post-qc"
readonly spark_dir="data/tmp/spark_dir"
readonly out_dir="data/mt/test"

readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_phased="${in_dir_phased}/ukb_wes_phased_non_singleton_chr${chr}-24xlong.qc-v4.2.2.vcf.gz"
readonly in_unphased="${in_dir_unphased}/ukb_wes_200k_filtered_chr${chr}.mt"

readonly out_prefix="${out_dir}/ukb_wes_200k_combined_chr${chr}"
readonly out="${out_prefix}.mt"

readonly hail_script="scripts/03_combine_mt.py"

mkdir -p ${out_dir}
if [ ! -f "${out}/_SUCCESS" ]; then
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy  
  warn_spark_dir_size
  python3 "${hail_script}" \
     --chrom ${chr} \
     --input_phased_path ${in_phased}\
     --input_unphased_path ${in_unphased} \
     --input_phased_type "vcf" \
     --input_unphased_type "mt" \
     --out_prefix ${out_prefix} \
     --out_type "mt"
  log_runtime ${SECONDS}
else
  print_update "file ${out} already exists. Skipping!"
fi



print_update "Finished running HAIL for chr${chr}" "${SECONDS}"





