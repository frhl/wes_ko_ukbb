#!/usr/bin/env bash

#SBATCH -A lindgren.prj
#SBATCH -J apply_int_to_cts
#SBATCH -D /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/apply_int_to_cts.log
#SBATCH --error=logs/apply_int_to_cts.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly r_script="scripts/phenotypes/00_apply_int_to_cts.R"
readonly in_dir="data/phenotypes"
readonly in_file="${in_dir}/curated_covar_phenotypes_cts_500k.tsv"

readonly out_dir="data/phenotypes"
readonly out_file="${out_dir}/curated_covar_phenotypes_cts_int_500k.txt.gz"

readonly phenotypes_cts=$(cat "${in_dir}/filtered_phenotypes_cts_manual.tsv" | tr "\n" "," | sed 's/\(.*\),/\1 /' )

set_up_rpy
Rscript ${r_script} \
  --in_file ${in_file} \
  --out_file ${out_file} \
  --phenotypes_list ${phenotypes_cts}






