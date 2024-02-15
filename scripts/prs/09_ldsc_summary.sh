#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=ldsc_summary
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/ldsc_summary.log
#SBATCH --error=logs/ldsc_summary.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --requeue
#SBATCH --begin=now+10hour

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/prs/09_ldsc_summary.R"

readonly ldsc_dir="data/prs/ldsc"
readonly pheno_dir="data/phenotypes"
readonly phenotypes="${pheno_dir}/dec22_phenotypes_binary_200k_header.tsv"

readonly out_dir="data/prs/validation"
readonly out_prefix="${out_dir}/ldsc_summary"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
    --in_dir "${ldsc_dir}" \
    --phenotypes "${phenotypes}" \
    --out_prefix "${out_prefix}" \
    --ldsc_n_eff_cutoff 5000 \
    --ldsc_p_cutoff 0.05


