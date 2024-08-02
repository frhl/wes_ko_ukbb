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

readonly rscript="scripts/saige/04_spa_combine.R"

# phenotypes we are including
readonly header_dir="data/phenotypes"
readonly header_file="${header_dir}/dec22_phenotypes_binary_200k_header.tsv"

# calc bonferroni cutoff
readonly genes_tested=952
readonly phenos_tested=311
readonly bonf_p_cutoff="$(python -c "print(0.05/(${genes_tested}*${phenos_tested}))")"
readonly nom_p_cutoff="$(python -c "print(0.05/(${genes_tested}))")"

# calc bonferonni cutoff for additive encoding
readonly add_genes_tested=16363
readonly add_phenos_tested=311
readonly add_bonf_p_cutoff="$(python -c "print(0.05/(${add_genes_tested}*${add_phenos_tested}))")"
readonly add_nom_p_cutoff="$(python -c "print(0.05/(${add_genes_tested}))")"


# define cutoffs
readonly N_ko_case_cutoff="0"
readonly N_ko_cutoff="5"

readonly out_dir="data/post_hoc/results/2024"
mkdir -p ${out_dir}

get_table() {
  local out_prefix="${out_dir}/${1}"
  local cond=${2}
  local prs=${3}
  local p=${4}
  if [ "${cond}" == "add_encoding" ]; then
    local gt_encoding="012"
  else 
    local gt_encoding="001"
  fi
  Rscript "${rscript}" \
    --path_header ${header_file} \
    --p_cutoff ${p} \
    --N_ko_case_cutoff ${N_ko_case_cutoff} \
    --N_ko_cutoff ${N_ko_cutoff} \
    --cond ${cond} \
    --prs ${prs} \
    --gt_encoding "${gt_encoding}" \
    --out_prefix "${out_prefix}"
}

set_up_rpy

# additive encoding
get_table "176k_sig_saige_all_prs_excl_N${N_ko_cutoff}_add_encoding" "add_encoding" "exclude" "1"
get_table "176k_sig_saige_all_prs_only_N${N_ko_cutoff}_add_encoding" "add_encoding" "only" "1"
get_table "176k_sig_saige_all_prs_pref_N${N_ko_cutoff}_add_encoding" "add_encoding" "prefer" "1"

#get_table "176k_sig_saige_sig_prs_excl_N${N_ko_cutoff}_add_encoding" "add_encoding" "exclude" "${add_nom_p_cutoff}"
get_table "176k_sig_saige_sig_prs_only_N${N_ko_cutoff}_add_encoding" "add_encoding" "only" "${add_nom_p_cutoff}"
get_table "176k_sig_saige_sig_prs_pref_N${N_ko_cutoff}_add_encoding" "add_encoding" "prefer" "${add_nom_p_cutoff}"

# Subsetting by P-value
get_table "176k_sig_saige_sig_prs_excl_N${N_ko_cutoff}" "none" "exclude" "${nom_p_cutoff}"
get_table "176k_sig_saige_sig_prs_only_N${N_ko_cutoff}" "none" "only" "${nom_p_cutoff}"
get_table "176k_sig_saige_sig_prs_pref_N${N_ko_cutoff}" "none" "prefer" "${nom_p_cutoff}"

# no subsetting by P-value
get_table "176k_sig_saige_all_prs_excl_N${N_ko_cutoff}" "none" "exclude" "1"
get_table "176k_sig_saige_all_prs_only_N${N_ko_cutoff}" "none" "only" "1"
get_table "176k_sig_saige_all_prs_pref_N${N_ko_cutoff}" "none" "prefer" "1"

# chet only
get_table "176k_sig_saige_all_prs_excl_N${N_ko_cutoff}_chet_only" "chet_only" "exclude" "1"
get_table "176k_sig_saige_all_prs_pref_N${N_ko_cutoff}_chet_only" "chet_only" "prefer" "1"

# hom only
get_table "176k_sig_saige_all_prs_excl_N${N_ko_cutoff}_hom_only" "hom_only" "exclude" "1"
get_table "176k_sig_saige_all_prs_pref_N${N_ko_cutoff}_hom_only" "hom_only" "prefer" "1"

# add encoding no pp cutoff
get_table "176k_sig_saige_all_prs_excl_N${N_ko_cutoff}_add_encoding_no_pp_cutoff" "add_encoding_no_pp_cutoff" "exclude" "1"
get_table "176k_sig_saige_all_prs_pref_N${N_ko_cutoff}_add_encoding_no_pp_cutoff" "add_encoding_no_pp_cutoff" "prefer" "1"

# no singletons
get_table "176k_sig_saige_all_prs_excl_N${N_ko_cutoff}_exclude_singletons" "exclude_singletons" "exclude" "1"
get_table "176k_sig_saige_all_prs_pref_N${N_ko_cutoff}_exclude_singletons" "exclude_singletons" "prefer" "1"

# only singletons (deprecated - not enough psuedo variants)
#get_table "176k_sig_saige_all_prs_excl_N${N_ko_cutoff}_only_singletons" "only_singletons" "exclude" "1"
#get_table "176k_sig_saige_all_prs_pref_N${N_ko_cutoff}_only_singletons" "only_singletons" "prefer" "1"








