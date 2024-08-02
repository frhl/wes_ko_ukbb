#!/usr/bin/env bash
#
#$ -N unphased_hardcalls
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/create_hardcalls_unphased.log
#$ -e logs/create_hardcalls_unphased.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 8
#$ -q short.qc@@short.hge
#$ -t 1-24

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories

readonly in_dir="data/unphased/post-qc"
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/hardcalls"

# hail script
readonly hail_script="scripts/QC/01_create_hardcalls.py"

# input path
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr${chr}.mt"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_chr${chr}"

mkdir -p ${out_dir}
if [ ! -f "${out}/_SUCCESS" ]; then
  set_up_hail
  set_up_vep
  set_up_pythonpath  
  python3 "${hail_script}" \
     --input_path ${in_file}\
     --input_type "mt" \
     --out_prefix ${out_prefix}
else
  print_update "file ${out} already exists. Skipping!"
fi


