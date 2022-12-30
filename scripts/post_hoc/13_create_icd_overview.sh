#!/usr/bin/env bash
#
#$ -N create_icd_overview
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/create_icd_overview.log
#$ -e logs/create_icd_overview.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/13_create_icd_overview.R"

readonly out_dir="data/knockouts/tables"
readonly out_prefix="${out_dir}/phenotype_icd_chapter"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
 --out_prefix "${out_prefix}" \



