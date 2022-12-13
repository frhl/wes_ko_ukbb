#!/usr/bin/env bash
#
#$ -N count_case_control_ko
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/count_case_control_ko.log
#$ -e logs/count_case_control_ko.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/post_hoc/02_count_case_ctrl_ko.R"

readonly pheno_dir="data/phenotypes"
readonly phenos="${pheno_dir}/spiros_brava_phenotypes_binary_200k.tsv.gz"
readonly ko_dir="data/knockouts/alt/pp90/combined"
readonly out_dir="data/post_hoc/results"

mkdir -p ${out_dir}

set_up_rpy

count_ko_case_ctrl() {
  local annotation="${1}"
  local ko_file="${ko_dir}/ukb_eur_wes_200k_chrCHR_${annotation}_all.tsv.gz"
  local out_prefix="${out_dir}/ukb_eur_wes_200k_case_ctrl_${annotation}"
  Rscript ${rscript} \
    --phenotypes ${phenos} \
    --ko_file ${ko_file} \
    --out_prefix ${out_prefix}
}


count_ko_case_ctrl "pLoF"
count_ko_case_ctrl "damaging_missense"
count_ko_case_ctrl "pLoF_damaging_missense"





