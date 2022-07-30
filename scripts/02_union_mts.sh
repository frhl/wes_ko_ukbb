#!/usr/bin/env bash
#
#$ -N union_mts
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/union_mts.log
#$ -e logs/union_mts.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 6
#$ -q long.qc@@long.hga
#$ -t 23

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark_dir"

# input files
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_dir_phased="data/phased/wes_union_calls/ligated"
readonly in_phased="${in_dir_phased}/ukb_eur_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly in_dir_unphased="data/unphased/wes/post-qc"
readonly in_unphased="${in_dir_unphased}/ukb_wes_200k_filtered_chr${chr}.mt"

# output files
readonly out_dir="data/mt/union"
readonly out_prefix="${out_dir}/ukb_eur_wes_union_calls_200k_chr${chr}"
readonly out="${out_prefix}.mt"

# hail script
readonly hail_script="scripts/02_union_mts.py"

mkdir -p ${out_dir}
if [ ! -f "${out}/_SUCCESS" ]; then
  if [ -f "${in_phased}" ]; then
    SECONDS=0
    set_up_hail
    set_up_pythonpath_legacy  
    python3 "${hail_script}" \
       --chrom ${chr} \
       --input_phased ${in_phased}\
       --input_unphased ${in_unphased} \
       --input_phased_type "vcf" \
       --input_unphased_type "mt" \
       --out_prefix ${out_prefix} \
       --out_type "mt" \
       && print_update "Finished merging phased and unphased data for chr${chr}" ${SECONDS} \
       || raise_error "Merging phased and unphased data for chr${chr} failed"
    log_runtime ${SECONDS}
  else
    >&2 echo "${in_phased} (input) does not exists. Exiting.."
  fi
else
  print_update "file ${out} already exists. Skipping!"
fi








