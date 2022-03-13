#!/usr/bin/env bash
#
#$ -N aggr_min_p
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/aggr_min_p.log
#$ -e logs/aggr_min_p.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1
#$ -V

#set -o errexit
#set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/permute/03_aggr_min_p.R"

readonly spa_cts_dir="data/saige/output/combined/cts/step2"
readonly spa_bin_dir="data/saige/output/combined/binary/step2"
readonly tsv_path="data/knockouts/alt_filtered/ukb_eur_wes_200k_chrCHR_maf0to5e-2_pLoF_damaging_missense.tsv.gz"
readonly out_dir="data/permute/overview"
readonly out_prefix="${out_dir}/overview"

mkdir -p ${out_dir}

set +eu
set_up_rpy
set -eu

Rscript ${rscript} \
  --tsv_path ${tsv_path} \
  --spa_cts_dir ${spa_cts_dir} \
  --spa_bin_dir ${spa_bin_dir} \
  --out_prefix ${out_prefix}



