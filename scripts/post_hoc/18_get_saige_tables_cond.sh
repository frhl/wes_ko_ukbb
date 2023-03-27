#!/usr/bin/env bash
#
#$ -N get_saige_tables_cond
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/get_saige_tables_cond.log
#$ -e logs/get_saige_tables_cond.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

#SBATCH --account=lindgren.prj
#SBATCH --job-name=get_saige_tables_cond
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/get_saige_tables_cond.log
#SBATCH --error=logs/get_saige_tables_cond.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/post_hoc/18_get_saige_tables_cond.R"

# phenotypes we are including
readonly header_dir="data/phenotypes"
readonly header_file="${header_dir}/dec22_phenotypes_binary_200k_header.tsv"

# for some reason there is a "bug" in saige, not giving us
# some of the columns when conditioniong on markers, 
# thus we need to map them manually
readonly ref_dir="data/post_hoc/results"
readonly ref_file="${ref_dir}/176k_sig_saige_sig_prs_pref.txt.gz"

# get merged hits that unprocessed
readonly merged_dir="data/conditional/combined/saige"
readonly merged_hits="${merged_dir}/176k_merged_hits_post_cond.txt.gz"

# paramters for sutff to include
readonly N_ko_case_cutoff="2"
readonly N_ko_cutoff="5"

readonly out_dir="data/post_hoc/results"
readonly out_prefix="${out_dir}/176k_sig_saige_cond_sig_pref_prs_combined"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --ref_file "${ref_file}" \
  --merged_hits "${merged_hits}" \
  --N_ko_case_cutoff ${N_ko_case_cutoff} \
  --N_ko_cutoff ${N_ko_cutoff} \
  --out_prefix "${out_prefix}"

