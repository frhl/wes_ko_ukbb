#!/usr/bin/env bash
#
#$ -N phased_hardcalls
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/create_hardcalls_phased.log
#$ -e logs/create_hardcalls_phased.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 8
#$ -q short.qc@@short.hge
#$ -t 21

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir="/well/lindgren/UKBIOBANK/nbaya/wes_200k/phase_ukb_wes/data/phased/non_singleton"
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/hardcalls/phased"

# hail script
readonly hail_script="scripts/QC/01_create_hardcalls.py"

# input path
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_file="${in_dir}/ukb_wes_phased_non_singleton_chr${chr}-24xlong.qc-v4.2.2.vcf.gz"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_chr${chr}"

mkdir -p ${out_dir}
if [ ! -f "${out}/_SUCCESS" ]; then
  set_up_hail
  set_up_vep
  set_up_pythonpath  
  python3 "${hail_script}" \
     --input_path ${in_file}\
     --input_type "vcf" \
     --out_prefix ${out_prefix}
else
  print_update "file ${out} already exists. Skipping!"
fi


