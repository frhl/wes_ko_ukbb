#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=get_attrition
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/get_attrition.log
#SBATCH --error=logs/get_attrition.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/saige/06_get_attrition.R"

#readonly dir="data/post_hoc/results"
readonly dir="data/post_hoc/results/2024"
readonly excl_prs="${dir}/176k_sig_saige_all_prs_excl_N5.txt.gz"
readonly prefer_prs="${dir}/176k_sig_saige_all_prs_pref_N5.txt.gz"
readonly cond_full="${dir}/176k_sig_saige_cond_all_pref_prs_combined.txt.gz"
readonly common_table="${dir}/176k_sig_saige_common_cond_sig_pref_prs_combined.txt.gz"

readonly additive="${dir}/176k_sig_saige_all_prs_pref_N5_add_encoding.txt.gz"
readonly cond_additive="${dir}/176k_sig_saige_dominance_cond_sig_pref_prs_combined.txt.gz"

readonly chet_only="${dir}/176k_sig_saige_all_prs_pref_N5_chet_only.txt.gz"
readonly hom_only="${dir}/176k_sig_saige_all_prs_pref_N5_hom_only.txt.gz"

readonly co_table_dir="data/knockouts/alt/pp90//co_occurence3"
readonly co_table="${co_table_dir}/co_occurence_collapsed_pLoF_damaging_missense.wide.txt.gz"

#readonly out_dir="data/post_hoc/results"
readonly out_dir="data/post_hoc/results/2024"
readonly out_prefix="${out_dir}/attrition"

set_up_rpy
Rscript "${rscript}" \
  --excl_prs "${excl_prs}" \
  --prefer_prs "${prefer_prs}" \
  --cond_full "${cond_full}" \
  --additive "${additive}" \
  --common "${common_table}" \
  --cond_additive "${cond_additive}" \
  --chet_only "${chet_only}" \
  --hom_only "${hom_only}" \
  --co_table "${co_table}" \
  --out_prefix "${out_prefix}"


