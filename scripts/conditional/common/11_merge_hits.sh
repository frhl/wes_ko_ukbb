#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=merge_hits
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/merge_hits.log
#SBATCH --error=logs/merge_hits.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly rscript="scripts/conditional/common/11_merge_hits.R"

readonly sig_genes_dir="data/conditional/combined/sig_genes/"
readonly sig_genes="${sig_genes_dir}/sig_genes_after_sig_prs_176k.txt.gz"
readonly in_dir="data/saige/output/binary/step2_common_new/min_mac4"
readonly out_dir="data/conditional/common/saige"
readonly out_prefix="${out_dir}/176k_merged_hits_post_common_cond"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --target_dir ${in_dir} \
  --path_sig_genes ${sig_genes} \
  --out_prefix ${out_prefix}



