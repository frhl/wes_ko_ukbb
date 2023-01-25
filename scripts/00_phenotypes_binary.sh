#!/usr/bin/env bash

#SBATCH -A lindgren.prj
#SBATCH -J phenotypes_binary
#SBATCH -D /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/phenotypes_binary.log
#SBATCH --error=logs/phenotypes_binary.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 4
#
#
#$ -N phenotypes_binary
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/phenotypes_binary.log
#$ -e logs/phenotypes_binary.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q short.qc
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly hail_script="scripts/00_phenotypes.py"
readonly r_script="scripts/00_phenotypes.R"
readonly spark_dir="data/tmp/spark"

readonly in_dir="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes"
readonly covar_dir="data/phenotypes"
readonly out_dir="data/phenotypes"

readonly in_bin="${in_dir}/curated_phenotypes_binary.tsv"
readonly in_cts="${in_dir}/curated_phenotypes_cts.tsv"
readonly tmp_bin="${out_dir}/dec22_phenotypes_binary.tsv.gz"
readonly out_bin_500k="${out_dir}/dec22_phenotypes_binary_500k"
readonly out_bin_200k="${out_dir}/new_dec22_phenotypes_binary_200k"

readonly vcf_sample_dir="data/phased/wes_union_calls/200k/shapeit5/ligated"
readonly vcf_sample="${vcf_sample_dir}/ukb_wes_union_calls_200k_chr22.vcf.bgz"

readonly final_sample_list="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list"
readonly final_ko_sample_list="data/phenotypes/samples/ukb_wes_ko.imputed.qc.samples"

readonly path_covars="${covar_dir}/covars1.csv"
readonly covariates="$(cat ${path_covars})"
mkdir -p ${out_dir}
mkdir -p ${spark_dir}

# Pre-processing of phenotypes
if [ ! -f ${tmp_bin} ]; then
  >&2 echo "Generating tmp files.."
  set_up_rpy
  Rscript ${r_script} \
    --input_path ${in_bin} \
    --covariates ${covariates} \
    --qc_samples ${final_sample_list} \
    --input_path_cts_to_bin ${in_cts} \
    --case_count_cutoff "50" \
    --include_spiros \
    --include_brava \
    --include_lindgren \
    --out_path ${tmp_bin}
fi

# set up python
set +eu
module purge
set_up_hail
set_up_pythonpath_legacy
set -eu

# get 500k WES samples
if [ ! -f ${out_bin_500k} ]; then
  >&2 echo "Genertating 500K.."
  python3 "${hail_script}" \
     --input_path "${tmp_bin}" \
     --export_header \
     --count_case_control \
     --out_prefix "${out_bin_500k}"
fi


if [ ! -f "${out_bin_500k}.tsv.gz" ]; then
  gzip "${out_bin_500k}.tsv"
fi

# Get 200k WES samples
if [ ! -f ${out_bin_200k} ]; then
  >&2 echo "Generating 200K.."
  python3 "${hail_script}" \
     --input_path "${tmp_bin}" \
     --extract_samples "${final_sample_list}" \
     --extract_phased_samples "${final_ko_sample_list}" \
     --export_header \
     --count_case_control \
     --out_prefix "${out_bin_200k}"
fi

if [ ! -f "${out_bin_200k}.tsv.gz" ]; then
  gzip "${out_bin_200k}.tsv"
fi








