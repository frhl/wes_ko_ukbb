#!/usr/bin/env bash

#SBATCH -A lindgren.prj
#SBATCH -J phenotypes_binary
#SBATCH -D /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/phenotypes_binary.log
#SBATCH --error=logs/phenotypes_binary.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --requeue

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly hail_script="scripts/00_phenotypes.py"
readonly r_script="scripts/00_phenotypes.R"
readonly spark_dir="data/tmp/spark"

readonly in_dir="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes"
readonly covar_dir="data/phenotypes"
readonly out_dir="data/phenotypes"

readonly in_bin="${in_dir}/curated_phenotypes_binary.tsv"
readonly tmp_bin="${out_dir}/spiros_phenotypes_binary.tsv.gz"
readonly out_bin_500k="${out_dir}/spiros_phenotypes_binary_500k"
readonly out_bin_200k="${out_dir}/spiros_phenotypes_binary_200k"

readonly final_sample_list="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list"

readonly path_covars="${covar_dir}/covars1.csv"
readonly covariates="$(cat ${path_covars})"

mkdir -p ${out_dir}

# Pre-processing of phenotypes
set_up_rpy
Rscript ${r_script} \
  --input_path ${in_bin} \
  --covariates ${covariates} \
  --qc_samples ${final_sample_list} \
  --include_spiros \
  --out_path ${tmp_bin}


# set up python
set +eu
module purge
set_up_hail
set_up_pythonpath_legacy
set -eu

# get 500k IMP samples
python3 "${hail_script}" \
     --input_path "${tmp_bin}" \
     --export_header \
     --count_case_control \
     --out_prefix "${out_bin_500k}"


# Get 200k WES samples
python3 "${hail_script}" \
     --input_path "${tmp_bin}" \
     --extract_samples "${final_sample_list}" \
     --export_header \
     --count_case_control \
     --out_prefix "${out_bin_200k}"






