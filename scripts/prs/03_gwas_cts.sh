#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=gwas_cts
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/gwas_cts.log
#SBATCH --error=logs/gwas_cts.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --array=1

set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly curwd=$(pwd)
readonly in_dir="data/prs/hapmap/ukb_500k/fitting"
readonly pheno_dir="data/phenotypes"

readonly bash_script="scripts/prs/_gwas_cts.sh"
readonly hail_script="scripts/prs/03_gwas_cts.py"
readonly merge_script="scripts/prs/_gwas_merge.sh"

readonly covar_file="${pheno_dir}/covars1.csv"
readonly covariates=$(cat ${covar_file})
#readonly covariates="age"

readonly input_type="plink"
readonly input_path="${in_dir}/ukb_hapmap_500k_eur_chrCHR"

readonly file_cts="${pheno_dir}/curated_covar_phenotypes_cts_int_500k.txt.gz"
readonly pheno_list_cts="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
readonly phenotype_cts=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${pheno_list_cts})

submit_gwas_job()
{
  local out_dir="${1}"
  local phenotype="${2}"
  local pheno_file="${3}"
  local slurm_tasks="${tasks}"
  local out_prefix="${out_dir}/ukb_hapmap_500k_eur_${phenotype}"
  local prefix="${out_prefix}_chrCHR"
  mkdir -p ${out_dir}
  if [ ! -f "${out_prefix}.txt.gz" ]; then
    if [ ! -z "${phenotype}" ] && [ ! -z "${covariates}" ]; then
      local gwas_jid=$(sbatch \
        --account="lindgren.prj" \
        --job-name="_gwas_cts_${phenotype}" \
        --output="logs/_gwas_cts.log" \
        --error="logs/_gwas_cts.errors.log" \
        --array=${slurm_tasks} \
        --chdir="$(pwd)" \
        --partition="short" \
        --cpus-per-task=3 \
        --parsable \
        "${bash_script}" \
        "${hail_script}" \
        "${input_path}" \
        "${input_type}" \
        "${pheno_file}" \
        "${phenotype}" \
        "${covariates}" \
        "${prefix}" )
      submit_merge_job "${out_dir}" "${prefix}" "${gwas_jid}"
    else
      >&2 echo "Required arguments are missing. Exiting.."
    fi
  else
    >&2 echo "${out_prefix}.txt.gz already exists. Skipping.."
  fi
}

submit_merge_job()
{
  local out_dir="${1}"
  local prefix="${2}"
  local gwas_jid="${3}"
  sbatch \
    --account="lindgren.prj" \
    --job-name="merge_cts_${phenotype}" \
    --output="${out_dir}/merge_${phenotype}.log" \
    --error="${out_dir}/merge_${phenotype}.errors.log" \
    --chdir="${curwd}" \
    --partition="short" \
    --cpus-per-task=1 \
    --dependency="afterok:${gwas_jid}" \
    --export=ALL \
    "${merge_script}" "${prefix}" "${out_dir}" "${out_prefix}.txt.gz"
}

readonly tasks="1-22"
submit_gwas_job "data/prs/sumstat/cts" "${phenotype_cts}" "${file_cts}"

