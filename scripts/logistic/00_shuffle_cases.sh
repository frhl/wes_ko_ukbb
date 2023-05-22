#!/usr/bin/env bash

#SBATCH -A lindgren.prj
#SBATCH -J shuffle_cases
#SBATCH -D /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/shuffle_cases.log
#SBATCH --error=logs/shuffle_cases.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly rscript="scripts/logistic/00_shuffle_cases.R"

readonly in_dir="data/phenotypes"
readonly in_file="${in_dir}/dec22_phenotypes_binary_200k.tsv.gz"
readonly seed=1

readonly out_dir="data/phenotypes/shuffle_cases"
readonly out_file="${out_dir}/dec22_phenotypes_binary_200k_shuffled_seed${seed}.txt.gz"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --input_path ${in_file} \
  --outfile ${out_file} \
  --seed ${seed}




