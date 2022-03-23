#!/usr/bin/env bash
#
#$ -N gwas
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/gwas.log
#$ -e logs/gwas.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 10-12
#$ -tc 2
#$ -V


set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly in_dir="data/prs/hapmap/ukb_500k"
readonly pheno_dir="data/phenotypes"

readonly bash_script="scripts/prs/_gwas.sh"
readonly hail_script="scripts/prs/03_gwas.py"
readonly merge_script="scripts/prs/_gwas_merge.sh"

readonly covar_file="${pheno_dir}/covars1.csv"
readonly covariates=$( cat ${covar_file} )

readonly input_type="plink"
readonly input_path="${in_dir}/ukb_hapmap_500k_eur_chrCHR"

readonly index=${SGE_TASK_ID}

readonly file_cts="${pheno_dir}/filtered_phenotypes_cts.txt" 
readonly pheno_list_cts="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
readonly phenotype_cts=$( sed "${index}q;d" ${pheno_list_cts} )

readonly file_binary="${pheno_dir}/filtered_phenotypes_binary.txt" 
readonly pheno_list_binary="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
readonly phenotype_binary=$( sed "${index}q;d" ${pheno_list_binary} )

readonly min_cases=100

submit_gwas_job()
{
  local out_dir="${1}"
  local phenotype="${2}"
  local pheno_file="${3}"
  local out_prefix="${out_dir}/ukb_hapmap_500k_eur_${phenotype}"
  local prefix="${out_prefix}_chrCHR"
  mkdir -p ${out_dir}
  set -x
  qsub -N "_${phenotype}_sumstat" \
    -t 1-22 \
    -tc 11 \
    -q short.qc@@short.hge \
    -pe shmem 1 \
    "${bash_script}" \
    "${hail_script}" \
    "${input_path}" \
    "${input_type}" \
    "${pheno_file}" \
    "${phenotype}" \
    "${covariates}" \
    "${min_cases}" \
    "${prefix}"
  set +x
  submit_merge_job
}

submit_merge_job()
{
  set -x
  qsub -N "_mrg_${phenotype}" \
    -q short.qc@@short.hge \
    -pe shmem 1 \
    -hold_jid "_${phenotype}_sumstat" \
    "${merge_script}" \
    "${prefix}" \
    "${out_dir}" \
    "${out_prefix}.txt.gz"
  set +x

}


#submit_gwas_job "data/prs/sumstat/binary" "${phenotype_binary}" "${file_binary}"
submit_gwas_job "data/prs/sumstat/cts" "${phenotype_cts}" "${file_cts}"



