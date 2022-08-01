#!/usr/bin/env bash
#
#$ -N calc_auc
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/calc_auc.log
#$ -e logs/calc_auc.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qf
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/prs/08_calc_auc.R"

readonly pgs_dir="data/prs/scores"
readonly out_dir="data/prs/validation"
readonly pheno_dir="data/phenotypes"

readonly phenotypes="${pheno_dir}/curated_covar_phenotypes_binary.tsv.gz"
readonly out_prefix="${out_dir}/pgs_auc_summary"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
    --directory "${pgs_dir}" \
    --phenotype_bin "${phenotypes}" \
    --out_prefix "${out_prefix}"


