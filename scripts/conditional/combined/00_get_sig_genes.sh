#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=get_sig_genes
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/get_sig_genes.log
#SBATCH --error=logs/get_sig_genes.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1

#$ -N get_sig_genes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/get_sig_genes.log
#$ -e logs/get_sig_genes.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qa
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly rscript="scripts/conditional/combined/00_get_sig_genes.R"

readonly cond_step="none" 
readonly p_cutoff="5e-7"
readonly prs="prefer"

readonly out_dir="data/conditional/combined/sig_genes"
readonly out_prefix="${out_dir}/sig_genes_after_prs_165k"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --cond_step ${cond_step} \
  --prs ${prs} \
  --p_cutoff ${p_cutoff} \
  --out_prefix ${out_prefix}



