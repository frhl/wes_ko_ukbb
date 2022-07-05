#!/usr/bin/env bash
#
#$ -N calc_cor
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/calc_cor.log
#$ -e logs/calc_cor.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qf
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/prs/09_calc_cor.R"

readonly pgs_dir="data/prs/scores"
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


