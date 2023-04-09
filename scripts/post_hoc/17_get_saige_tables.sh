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

#SBATCH --account=lindgren.prj
#SBATCH --job-name=get_saige_tables
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/get_saige_tables.log
#SBATCH --error=logs/get_saige_tables.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/post_hoc/17_get_saige_tables.R"

# phenotypes we are including
readonly header_dir="data/phenotypes"
readonly header_file="${header_dir}/dec22_phenotypes_binary_200k_header.tsv"

# calc bonferroni cutoff
readonly genes_tested=958
readonly phenos_tested=311
readonly p_cutoff="$(python -c "print(0.05/(${genes_tested}*${phenos_tested}))")"

# define cutoffs
readonly N_ko_case_cutoff="2"
readonly N_ko_cutoff="5"

# path to PRS heritability estimates
readonly ldsc_h2_dir="data/prs/validation"
readonly ldsc_h2="${ldsc_h2_dir}/ldsc_summary.txt.gz"

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
    --path_ldsc_h2 ${ldsc_h2} \
    --N_ko_case_cutoff ${N_ko_case_cutoff} \
    --N_ko_cutoff ${N_ko_cutoff} \
    --cond ${cond} \
    --prs ${prs} \
    --out_prefix "${out_prefix}"
}

set_up_rpy
# Subsetting by P-value
get_table "176k_sig_saige_sig_prs_excl" "none" "exclude" "${p_cutoff}"
get_table "176k_sig_saige_sig_prs_only" "none" "only" "${p_cutoff}"
get_table "176k_sig_saige_sig_prs_pref" "none" "prefer" "${p_cutoff}"

# no subsetting by P-value
get_table "176k_sig_saige_all_prs_excl" "none" "exclude" "1"
get_table "176k_sig_saige_all_prs_only" "none" "only" "1"
get_table "176k_sig_saige_all_prs_pref" "none" "prefer" "1"



