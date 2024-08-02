#!/usr/bin/env bash
#
#$ -N sample_geneset_counts
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/sample_geneset_counts.log
#$ -e logs/sample_geneset_counts.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qc

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/06_sample_geneset_counts.R"

readonly out_dir="data/knockouts/tables"
readonly out_prefix="${out_dir}/sample_geneset_counts"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --out_prefix "${out_prefix}"



