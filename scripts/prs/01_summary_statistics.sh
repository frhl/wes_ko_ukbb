#!/usr/bin/env bash
#
#$ -N summary_statistics
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/summary_statistics.log
#$ -e logs/summary_statistics.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hga
#$ -t 1-22
#$ -V

# cts
# 1 - 33 contains non-residuals
# 34 - 103 contains residuals

module purge
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly in_dir="data/prs/hapmap"
readonly pheno_dir="data/phenotypes"
readonly out_dir="data/prs/sumstat"

readonly bash_script="scripts/prs/_summary_statistics.sh"
readonly hail_script="scripts/prs/01_summary_statistics.py"
readonly merge_script="scripts/prs/_summary_statistics_merge.sh"

readonly covar_file="${pheno_dir}/covars1.csv"
readonly covariates=$( cat ${covar_file} )

readonly input_type="plink"
readonly input_path="${in_dir}/ukb_hapmap_500k_eur_chrCHR"
readonly pheno_file="${pheno_dir}/curated_phenotypes.tsv" 

readonly pheno_list_cts="${pheno_dir}/curated_phenotypes_cts_header.tsv"
readonly phenotype_cts=$( cut -f${SGE_TASK_ID} ${pheno_list_cts} )
readonly pheno_list_binary="${pheno_dir}/curated_phenotypes_binary_header.tsv"
readonly phenotype_binary=$( cut -f${SGE_TASK_ID} ${pheno_list_binary} )

#readonly out_prefix="${out_dir}/ukb_hapmap_500k_eur_${phenotype}"
#readonly prefix="${out_prefix}_chrCHR"

submit_sumstat_job()
{
  mkdir -p ${out_dir}
  local phenotype="${1}"
  local out_prefix="${out_dir}/ukb_hapmap_500k_eur_${phenotype}"
  local prefix="${out_prefix}_chrCHR"
  set -x
  qsub -N "_${phenotype}_sumstat" \
    -t 1-22 \
    -q short.qc@@short.hge \
    -pe shmem 1 \
    "${bash_script}" \
    "${hail_script}" \
    "${input_path}" \
    "${input_type}" \
    "${pheno_file}" \
    "${phenotype}" \
    "${covariates}" \
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

submit_sumstat_job "${phenotype_binary}"



