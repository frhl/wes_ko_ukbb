#!/usr/bin/env bash
#
#$ -N ukb_500k_binary
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/ukb_500k_binary.log
#$ -e logs/ukb_500k_binary.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly r_script="scripts/prs/00_phenotypes_500k.R"

readonly in_dir="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes"
readonly covar_dir="data/phenotypes"
readonly out_dir="data/phenotypes/500k"

readonly in_bin="${in_dir}/curated_phenotypes_binary.tsv"

readonly out_bin="${out_dir}/500k_curated_covar_phenotypes_binary.tsv.gz"

readonly path_covars="${covar_dir}/covars1.csv"
readonly covariates="$(cat ${path_covars})"

mkdir -p ${out_dir}

# * add age2, age3 and sex-age covariates
# * add sex-stratified residiuals"
set_up_rpy
Rscript ${r_script} \
  --input_path ${in_bin} \
  --out_path ${tmp_bin}




