#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=get_attrition
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/get_attrition.log
#SBATCH --error=logs/get_attrition.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/13_get_attrition.R"

readonly header_dir="data/phenotypes"
readonly header_file="${header_dir}/dec22_phenotypes_binary_200k_header.tsv"

set_up_rpy
Rscript "${rscript}" \
  --ref_file "${ref_file}" \
  --merged_hits "${merged_hits}" \
  --N_ko_case_cutoff ${N_ko_case_cutoff} \
  --N_ko_cutoff ${N_ko_cutoff} \
  --out_prefix "${out_prefix}"


