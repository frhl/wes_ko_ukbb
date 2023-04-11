#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=spa_combine
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/spa_combine.log
#SBATCH --error=logs/spa_combine.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/12_spa_combine.R"

# phenotypes we are including
readonly header_dir="data/phenotypes"
readonly header_file="${header_dir}/dec22_phenotypes_binary_200k_header.tsv"

# calc bonferroni cutoff
readonly genes_tested=958
readonly phenos_tested=311
readonly p_cutoff="$(python -c "print(0.05/(${genes_tested}*${phenos_tested}))")"

# define cutoffs
readonly N_ko_case_cutoff="0"
readonly N_ko_cutoff="5"

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
# Subsetting by P-value
#get_table "176k_sig_saige_sig_prs_excl_wo_case_cutoff" "none" "exclude" "${p_cutoff}"
#get_table "176k_sig_saige_sig_prs_only_wo_case_cutoff" "none" "only" "${p_cutoff}"
#get_table "176k_sig_saige_sig_prs_pref_wo_case_cutoff" "none" "prefer" "${p_cutoff}"

# no subsetting by P-value
get_table "176k_sig_saige_all_prs_excl_wo_case_cutoff" "none" "exclude" "1"
get_table "176k_sig_saige_all_prs_only_wo_case_cutoff" "none" "only" "1"
get_table "176k_sig_saige_all_prs_pref_wo_case_cutoff" "none" "prefer" "1"



