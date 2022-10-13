#!/usr/bin/env bash
#
# @description Annotate main MatrixTables with VEP results
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=annotate
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/annotate.log
#SBATCH --error=logs/annotate.errors.log
#SBATCH --partition=long
#SBATCH --cpus-per-task 2
#SBATCH --array=1-22
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_dir="data/mt/annotated"
readonly input_prefix="${in_dir}/ukb_eur_wes_union_calls_200k_chr${chr}.mt"

mkdir -p ${out_dir}

python3 "${hail_script}" \
   --input_path ${input_prefix}\
   --input_type "mt" \






