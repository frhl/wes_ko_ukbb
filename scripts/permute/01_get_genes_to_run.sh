#!/usr/bin/env bash
#
# @description for each gene across all phenotypes figure out how many permutations are needed
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=get_genes_to_run
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/get_genes_to_run.log
#SBATCH --error=logs/get_genes_to_run.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/permute/01_get_genes_to_run.R"

readonly min_chet=2
readonly min_cis=2

readonly out_dir="data/permute/overview"
readonly out_prefix="${out_dir}/genes_to_run_${min_cis}cis_${min_chet}chets"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --out_prefix ${out_prefix} \
  --min_chet ${min_chet} \
  --min_cis ${min_cis}



