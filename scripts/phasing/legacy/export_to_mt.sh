#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=export_to_mt
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/export_to_mt.log
#SBATCH --error=logs/export_to_mt.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=1-22
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="scripts/phasing/export_to_mt.py"

readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} ) 
readonly in_dir="data/unphased/wes_union_calls"
readonly input_prefix="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly input_type="vcf"

readonly out_dir="data/unphased/wes_union_calls"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}"
readonly out_type="mt"


mkdir -p ${out_dir}

set_up_hail
set_up_pythonpath_legacy  
python3 "${hail_script}" \
   --in_file "${input_prefix}" \
   --in_type "${input_type}" \
   --out_prefix ${out_prefix} \
   --out_type "${out_type}" \


