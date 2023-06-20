#!/usr/bin/env bash
#
# Iterate over files in data/saige/output/binary/step2 and combine
# combine significant markers into a single file
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=get_sig_genes
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/get_sig_genes.log
#SBATCH --error=logs/get_sig_genes.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly rscript="scripts/conditional/combined/00_get_sig_genes.R"

readonly cond_step="none" 
readonly p_cutoff="5.2e-5" # nominal sig
readonly prs="prefer"

readonly out_dir="data/conditional/combined/sig_genes"
readonly out_prefix="${out_dir}/test_sig_genes_after_sig_prs_176k"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --cond_step ${cond_step} \
  --prs ${prs} \
  --p_cutoff ${p_cutoff} \
  --out_prefix ${out_prefix}



