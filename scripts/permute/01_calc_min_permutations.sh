#!/usr/bin/env bash
#
# @description for each gene across all phenotypes figure out how many permutations are needed
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=calc_min_permutations
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/calc_min_permutations.log
#SBATCH --error=logs/calc_min_permutations.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/permute/01_calc_min_permutations.R"

# set parameters
readonly min_mac=4
readonly p_cutoff="5e-7"

# directories and out paths
readonly basename_prefix="ukb_eur_wes_200k_maf0to5e-2"
readonly spa_cts_dir="data/saige/output/cts/step2_common_cond/min_mac${min_mac}"
readonly spa_bin_dir="data/saige/output/binary/step2_common_cond/min_mac${min_mac}"
readonly out_dir="data/permute/overview/min_mac${min_mac}/phased_only"
readonly out_prefix="${out_dir}/main"
mkdir -p ${out_dir}

# required to pin down which genes are ch knockout and which are hom alt knockouts
readonly tsv_path="data/knockouts/alt/ukb_eur_wes_200k_chrCHR_maf0to5e-2_pLoF_damaging_missense.tsv.gz"

set_up_rpy
Rscript ${rscript} \
  --tsv_path ${tsv_path} \
  --spa_cts_dir ${spa_cts_dir} \
  --spa_bin_dir ${spa_bin_dir} \
  --basename_prefix ${basename_prefix} \
  --out_prefix ${out_prefix} \
  --use_cond_p



