#!/usr/bin/env bash
#
# @description for each gene across all phenotypes figure out how many permutations are needed
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=calc_min_permutations
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/calc_min_permutations.log
#SBATCH --error=logs/calc_min_permutations.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/permute/01_calc_min_permutations.R"

# set parameters
readonly min_mac=4
readonly p_cutoff="5e-7"

# directories and out paths
readonly out_dir="data/permute/overview/min_mac${min_mac}/no_cond"
readonly out_prefix="${out_dir}/main"
mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --out_prefix ${out_prefix} \
  --cond "none" \
  --prs "prefer"
  #--use_cond_p



