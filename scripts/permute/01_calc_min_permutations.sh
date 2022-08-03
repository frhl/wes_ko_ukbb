#!/usr/bin/env bash
#
#$ -N calc_min_permutations
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/calc_min_permutations.log
#$ -e logs/calc_min_permutations.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/permute/01_calc_min_permutations.R"

# set parameters
readonly min_mac=4
readonly p_cutoff="5e-7"

# directories and out paths
readonly spa_cts_dir="data/saige/output/cts/step2_common_cond/min_mac${min_mac}"
readonly spa_bin_dir="data/saige/output/binary/step2_common_cond/min_mac${min_mac}"
readonly out_dir="data/permute/overview/min_mac${min_mac}"
readonly out_prefix="${out_dir}/main"
mkdir -p ${out_dir}

# required to pin down which genes are ch knockout and which are hom alt knockouts
readonly tsv_path="data/knockouts/alt/ukb_eur_wes_200k_chrCHR_maf0to5e-2_pLoF_damaging_missense.tsv.gz"

set_up_rpy
Rscript ${rscript} \
  --tsv_path ${tsv_path} \
  --spa_cts_dir ${spa_cts_dir} \
  --spa_bin_dir ${spa_bin_dir} \
  --out_prefix ${out_prefix} \
  --use_cond_p



