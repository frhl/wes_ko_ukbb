#!/usr/bin/env bash
#
#$ -N calc_min_permutations
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/calc_min_permutations.log
#$ -e logs/calc_min_permutations.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1
#$ -V

#set -o errexit
#set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/permute/01_calc_min_permutations.R"

readonly spa_cts_dir="data/saige/output/combined/cts/step2"
readonly spa_bin_dir="data/saige/output/combined/binary/step2"
readonly out_dir="data/permute/overview"
readonly out_prefix="${out_dir}/overview"

# required to pin down which genes are ch knockout and which are hom alt knockouts
readonly tsv_path="data/knockouts/alt_filtered/ukb_eur_wes_200k_chrCHR_maf0to5e-2_pLoF_damaging_missense.tsv.gz"

mkdir -p ${out_dir}

set_up_rpy

Rscript ${rscript} \
  --tsv_path ${tsv_path} \
  --spa_cts_dir ${spa_cts_dir} \
  --spa_bin_dir ${spa_bin_dir} \
  --out_prefix ${out_prefix}



