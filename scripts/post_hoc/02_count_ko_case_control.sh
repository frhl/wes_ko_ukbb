#!/usr/bin/env bash
#
#$ -N count_ko_cases_controls
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/count_ko_cases_controls.log
#$ -e logs/count_ko_cases_controls.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/permute/02_count_ko_cases_controls.R"

readonly pheno_dir="data/phenotypes"
readonly phenos="${pheno_dir}/filtered_phenotypes_binary.tsv"
readonly ko_dir="data/knockouts/alt"
readonly ko_file="${ko_dir}/ukb_eur_wes_200k_chrCHR_maf0to5e-2_damaging_missense.tsv.gz"
readonly out_dir="data/post_hoc/results"
readonly out_prefix="${out_dir}/count_ko_case_control"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --phenotypes ${phenos} \
  --ko_file ${ko_file} \
  --out_prefix ${out_prefix}








