#!/usr/bin/env bash
#
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=calc_cor
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/calc_cor.log
#SBATCH --error=logs/calc_cor.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/prs/09_calc_cor.R"

readonly pgs_dir="data/prs/scores_full"
readonly out_dir="data/prs/validation"
readonly pheno_dir="data/phenotypes"

readonly phenotypes="${pheno_dir}/curated_covar_phenotypes_cts.tsv.gz"
readonly out_prefix="${out_dir}/pgs_cor_summary"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
    --directory "${pgs_dir}" \
    --phenotype_cts "${phenotypes}" \
    --out_prefix "${out_prefix}"


