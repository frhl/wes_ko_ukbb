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

readonly rscript="scripts/post_hoc/04_create_icd_overview.R"

readonly in_dir="data/phenotypes"
readonly in_file="${in_dir}/dec22_phenotypes_binary_200k_header.tsv"
readonly desc_short="${in_dir}/phenotype_icd_chapter_desc_short.txt.gz"

readonly out_dir="data/phenotypes"
readonly out_prefix="${out_dir}/phenotype_icd_chapter"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --out_prefix "${out_prefix}" \
  --phenotypes_to_keep "${in_file}" \
  --path_icd_desc_short "${desc_short}"

