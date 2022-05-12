#!/usr/bin/env bash
#
#$ -N collect_permute
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/collect_permute.log
#$ -e logs/collect_permute.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/permute/05_collect_permute.R"

readonly target_dir="data/permute/permutations"
readonly out_dir="data/permute/results"
readonly out_prefix="${out_dir}/counts"


set_up_rpy
Rscript ${rscript} \
  --target_dir ${target_dir} \
  --out_prefix ${out_prefix}








