#!/usr/bin/env bash
#
# We only want to permute phenotypes in which there is a significant assocation
# and enough chets & multi-hit cis
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=get_phenos_to_run
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/get_phenos_to_run.log
#SBATCH --error=logs/get_phenos_to_run.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/permute/00_get_phenos_to_run.R"

readonly min_chet=2
readonly min_cis=2

readonly sig_hits_dir="data/post_hoc/results/N_ko10"
readonly sig_hits="${sig_hits_dir}/176k_sig_saige_sig_prs_pref_N5.txt.gz"

readonly out_dir="data/permute/overview/N_ko10"
readonly out_prefix="${out_dir}/phenotypes_with_${min_cis}cis_${min_chet}chets"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --min_chet ${min_chet} \
  --min_cis ${min_cis} \
  --path_sig_hits ${sig_hits} \
  --out_prefix ${out_prefix}




