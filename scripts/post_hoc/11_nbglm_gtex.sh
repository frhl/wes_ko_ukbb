#!/usr/bin/env bash
#
#$ -N nbglm_gtex
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/nbglm_gtex.log
#$ -e logs/nbglm_gtex.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/11_nbglm_gtex.R"

readonly in_dir="data/knockouts/tables"
readonly in_file="${in_dir}/combined_annotations_by_sample.counts.txt.gz"

readonly out_dir="data/knockouts/tables"
readonly out_prefix="${out_dir}/nbglm_gtex"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --out_prefix "${out_prefix}" \
  --consolidated "${in_file}"


