#!/usr/bin/env bash
#
#$ -N get_saige_tables
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/get_saige_tables.log
#$ -e logs/get_saige_tables.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/post_hoc/17_get_saige_tables.R"

# phenotypes we are including
readonly header_dir="data/phenotypes"
readonly header_file="${header_dir}/dec22_phenotypes_binary_200k_header.tsv"

# bonferroni p-value < 0.05 / ( 313 * 1895 ) # when N_ko>=2 (min_mac=4)
# bonferroni p-value < 0.05 / (313 * 1143) # when N_ko>=2 (min_mac=8)
readonly p_cutoff="1.3975888796648e-07"
readonly N_ko_case_cutoff="1"
readonly N_ko_cutoff="4"

readonly out_dir="data/post_hoc/results"
mkdir -p ${out_dir}

get_table() {
  local out_prefix="${out_dir}/${1}"
  local cond=${2}
  local prs=${3}
  local p=${4}
  Rscript "${rscript}" \
    --path_header ${header_file} \
    --p_cutoff ${p} \
    --N_ko_case_cutoff ${N_ko_case_cutoff} \
    --N_ko_cutoff ${N_ko_cutoff} \
    --cond ${cond} \
    --prs ${prs} \
    --out_prefix "${out_prefix}"
}

set_up_rpy
get_table "165k_saige_sig_merge_excluding_prs" "none" "exclude" "${p_cutoff}"
get_table "165k_saige_sig_merge_only_prs" "none" "only" "${p_cutoff}"
get_table "165k_saige_sig_merge_prefer_prs" "none" "prefer" "${p_cutoff}"

get_table "165k_saige_merge_excluding_prs" "none" "exclude" "1"
get_table "165k_saige_merge_only_prs" "none" "only" "1"
get_table "165k_saige_merge_prefer_prs" "none" "prefer" "1"



