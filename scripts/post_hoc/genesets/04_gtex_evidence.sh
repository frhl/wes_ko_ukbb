#!/usr/bin/env bash
#
#$ -N gtex_evidence
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/gtex_evidence.log
#$ -e logs/gtex_evidence.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/04_gtex_evidence.R"
readonly out_dir="derived/tables"
readonly out_prefix="${out_dir}/gtex_top10"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --out_prefix ${out_prefix}


