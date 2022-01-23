#!/usr/bin/env bash
#
#$ -N summary_statistics
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/summary_statistics.log
#$ -e logs/summary_statistics.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hga
#$ -t 34
#$ -V

# 1 - 33 contains non-residuals
# 34 - 103 contains residuals

module purge
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly pheno_dir="data/phenotypes"
readonly out_dir="data/prs/sumstat"

readonly bash_script="scripts/prs/_summary_statistics.sh"
readonly hail_script="scripts/prs/01_summary_statistics.py"

readonly covar_file="${pheno_dir}/covars1.csv"
readonly covariates=$( cat ${covar_file} )

readonly pheno_file="${pheno_dir}/curated_phenotypes.tsv" 
readonly pheno_list="${pheno_dir}/curated_phenotypes_cts_header.tsv"
readonly phenotype=$( cut -f${SGE_TASK_ID} ${pheno_list} )

readonly out_prefix="${out_dir}/ukb_imp_500k_${phenotype}"

readonly dataset="imp"
readonly liftover="yes"
readonly min_info="0.80"

submit_sumstat_job()
{
  mkdir -p ${out_dir}
  local prefix="${out_prefix}_chrCHR"
  set -x
  qsub -N "_${phenotype}_sumstat" \
    -t 21 \
    -q short.qc@@short.hga \
    -pe shmem 1 \
    "${bash_script}" \
    "${hail_script}" \
    "${dataset}" \
    "${pheno_file}" \
    "${phenotype}" \
    "${response}" \
    "${covariates}" \
    "${min_info}" \
    "${prefix}" \
    "${liftover}"
  set +x
}

submit_sumstat_job

