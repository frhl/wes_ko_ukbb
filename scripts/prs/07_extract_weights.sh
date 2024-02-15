#!/usr/bin/env bash
#
# Fit weights for PGS
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=extract_weights
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/extract_weights.log
#SBATCH --error=logs/extract_weights.errors.log
#SBATCH --open-mode=append
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-100
#
#$ -N extract_weights
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/extract_weights.log
#$ -e logs/extract_weights.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1-320
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/prs/07_extract_weights.R"

readonly ldsc_dir="data/prs/ldsc"
readonly weights_dir="data/prs/weights/auto"
readonly pheno_dir="data/phenotypes"
readonly out_dir="data/prs/weights"

readonly cluster=$( get_current_cluster )
readonly index=$( get_array_task_id )

readonly file_binary="${pheno_dir}/dec22_phenotypes_binary_500k.tsv.gz"
readonly pheno_list_binary="${pheno_dir}/dec22_phenotypes_binary_200k_header.tsv"
readonly phenotype_binary=$( sed "${index}q;d" ${pheno_list_binary} )

mkdir -p ${out_dir}

extract_weights() {
  set_up_ldpred2
  local phenotype=${1}
  local ldsc="${ldsc_dir}/ldsc_${phenotype}.rds"
  local path_betas="${weights_dir}/weights_${phenotype}.rda"
  local path_betas_map="${weights_dir}/weights_${phenotype}.txt.gz"
  local out_prefix="${out_dir}/weights.qc.${phenotype}"
  if [ ! -z ${phenotype} ]; then
     if [ ! -f "${out_prefix}.txt.gz" ]; then
          Rscript "${rscript}" \
            --ldsc "${ldsc}" \
            --path_betas "${path_betas}" \
            --path_betas_map "${path_betas_map}" \
            --out_prefix "${out_prefix}"
     else
        >&2 echo "${out_prefix}.txt.gz already exists. Skipping.."
     fi
  fi
}

extract_weights ${phenotype_binary}


