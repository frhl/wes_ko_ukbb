#!/usr/bin/env bash
#
# @description for each gene across all phenotypes figure out how many permutations are needed
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=get_min_p
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/get_min_p.log
#SBATCH --error=logs/get_min_p.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/permute/01_get_min_p.R"

# set parameters
readonly min_mac=4
readonly p_cutoff="5e-7"

# only permute genes from phentypes where we 
# have a chance of detecting effecs
readonly in_dir="data/permute/overview/min_mac4"
readonly in_file="${in_dir}/phenotypes_with_5cis_5chets.txt.gz"

# subset phenotypes
readonly phenotypes=$( zcat ${in_file} | cut -f1 | tail -n +2 | tr "\n" ",")

# directories and out paths
readonly out_dir="data/permute/overview/min_mac${min_mac}"
readonly out_prefix="${out_dir}/main"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --out_prefix ${out_prefix} \
  --phenotypes ${phenotypes} \
  --use_cond_p \
  --cond "none" \
  --prs "prefer"



