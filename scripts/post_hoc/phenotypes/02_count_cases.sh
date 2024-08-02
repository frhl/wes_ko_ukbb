#!/usr/bin/env bash
#
#$ -N count_cases
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/count_cases.log
#$ -e logs/count_cases.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/post_hoc/02_count_cases.R"

readonly pheno_dir="data/phenotypes"
readonly phenos_boolean="${pheno_dir}/dec22_phenotypes_binary_200k.tsv.gz"
readonly phenos_header="${pheno_dir}/dec22_phenotypes_binary_200k_header.tsv"
readonly phenos_tte="${pheno_dir}/tte_matrix_176k.txt.gz"

readonly sample_dir="data/phenotypes/samples"
readonly samples="${sample_dir}/ukb_wes_ko.qc.nfe.samples"

readonly out_dir="data/post_hoc/results"

mkdir -p ${out_dir}

set_up_rpy

count_ko_case_ctrl() {
  local phenos="${1}"
  local out_tag="${2}"
  local annotation="${3}"
  local out_prefix="${out_dir}/ukb_eur_wes_200k_case_ctrl_${annotation}_${out_tag}"
  Rscript ${rscript} \
    --phenotypes ${phenos} \
    --phenos_path ${phenos_header} \
    --annotation ${annotation} \
    --samples ${samples} \
    --out_prefix ${out_prefix}
}

count_ko_case_ctrl "${phenos_tte}" "tte" "pLoF"
#count_ko_case_ctrl "${phenos_boolean}" "bool" "pLoF"
#count_ko_case_ctrl "damaging_missense"
#count_ko_case_ctrl "pLoF_damaging_missense"
#count_ko_case_ctrl "other_missense"
#count_ko_case_ctrl "synonymous"





