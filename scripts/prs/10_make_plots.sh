#!/usr/bin/env bash
#
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=make_plots
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/make_plots.log
#SBATCH --error=logs/make_plots.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/prs/10_make_plots.R"


readonly pgs_dir="data/prs/validation"
readonly out_dir="data/prs/validation"
readonly pheno_dir="data/phenotypes"

readonly summary_ldsc="${pgs_dir}/ldsc_summary.txt.gz"
readonly summary_prs_cts="${pgs_dir}/pgs_cor_summary.txt.gz"
readonly summary_prs_bin="${pgs_dir}/pgs_auc_summary.txt.gz"

readonly bin_header="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
readonly cts_header="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"

readonly phenotypes="${pheno_dir}/curated_covar_phenotypes_cts.tsv.gz"
readonly out_prefix="${out_dir}/220712_pgs_results"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
    --out_prefix "${out_prefix}" \
    --summary_ldsc "${summary_ldsc}" \
    --summary_prs_cts "${summary_prs_cts}" \
    --summary_prs_bin "${summary_prs_bin}" \
    --bin_header "${bin_header}" \
    --cts_header "${cts_header}"



