#!/usr/bin/env bash
#
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=gwas
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/gwas.log
#SBATCH --error=logs/gwas.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=90
#SBATCH --requeue


set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly curwd=$(pwd)
readonly in_dir="data/prs/hapmap/ukb_500k/fitting"
readonly pheno_dir="data/phenotypes"

readonly bash_script="scripts/prs/_gwas.sh"
readonly hail_script="scripts/prs/03_gwas.py"
readonly merge_script="scripts/prs/_gwas_merge.sh"

readonly covar_file="${pheno_dir}/covars1.csv"
readonly covariates=$( cat ${covar_file} )

readonly input_type="plink"
readonly input_path="${in_dir}/ukb_hapmap_500k_eur_chrCHR"

readonly index=${SLURM_ARRAY_TASK_ID}

readonly file_cts="${pheno_dir}/curated_covar_phenotypes_cts.tsv.gz" 
readonly pheno_list_cts="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
readonly phenotype_cts=$( sed "${index}q;d" ${pheno_list_cts} )

readonly file_binary="${pheno_dir}/spiros_brava_phenotypes_binary_200k.tsv"
readonly pheno_list_binary="${pheno_dir}/spiros_brava_phenotypes_binary_200k_header.tsv"
readonly phenotype_binary=$( sed "${index}q;d" ${pheno_list_binary} )

# the job will fail if less than X cases
readonly min_cases=100

submit_gwas_job()
{
  local out_dir="${1}"
  local phenotype="${2}"
  local pheno_file="${3}"
  local out_prefix="${out_dir}/ukb_hapmap_500k_eur_${phenotype}"
  local prefix="${out_prefix}_chrCHR"
  mkdir -p ${out_dir}
  if [ ! -f "${out_prefix}.txt.gz" ]; then
    if [ ! -z ${phenotype} ]; then
      if [ ! -z ${covariates} ]; then
        local slurm_jname="_${phenotype}_gwas"
        local slurm_lname="logs/_gwas"
        local slurm_project="lindgren.prj"
        local slurm_tasks="${tasks}"
        local slurm_queue="short"
        local slurm_shmem="3"
        readonly gwas_jid=$( sbatch \
          --account="${slurm_project}" \
          --job-name="${slurm_jname}" \
          --output="${slurm_lname}.log" \
          --error="${slurm_lname}.errors.log" \
          --chdir="$(pwd)" \
          --partition="${slurm_queue}" \
          --cpus-per-task="${slurm_shmem}" \
          --array=${slurm_tasks} \
          --parsable \
          "${bash_script}" \
          "${hail_script}" \
          "${input_path}" \
          "${input_type}" \
          "${pheno_file}" \
          "${phenotype}" \
          "${covariates}" \
          "${min_cases}" \
          "${prefix}" )
        submit_merge_job
      else
      >&2 echo "Covariate argument is emtpy! Exiting.."
      fi
    else
      >&2 echo "Phenotype argument is empty. Exiting.."  
    fi
  else
    >&2 echo "${out_prefix}.tsv.gz already exists. Skipping.."
  fi
}

submit_merge_job()
{
  local slurm_jname="_mrg_${phenotype}"
  local slurm_lname="_gwas_merge"
  local slurm_project="lindgren.prj"
  local slurm_queue="short"
  local slurm_nslots="1"
  readonly merge_jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname}.log" \
    --error="${slurm_lname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --dependency="afterok:${gwas_jid}" \
    --open-mode="append" \
    --parsable \
    "${merge_script}" \
    "${prefix}" \
    "${out_dir}" \
    "${out_prefix}.txt.gz" )
}


readonly tasks="1-22"
submit_gwas_job "data/prs/sumstat/binary" "${phenotype_binary}" "${file_binary}"
#submit_gwas_job "data/prs/sumstat/cts" "${phenotype_cts}_int" "${file_cts}"



