#!/usr/bin/env bash
#
# extract genotypes in regions near genes that are significant in primary analysis
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=filter_genotypes
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/filter_genotypes.log
#SBATCH --error=logs/filter_genotypes.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-10

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly curwd=$(pwd)
readonly bash_script="scripts/conditional/common/_filter_genotypes.sh"

readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'

# padding upstream/downstream
readonly padding=1000000 # 1 Megabase
# minimum maf and imputation score extracted
readonly min_maf=0.01
readonly min_info=0.8
readonly min_mac=4

readonly in_dir="data/conditional/common/gene_positions/min_mac${min_mac}"
readonly out_dir="data/conditional/common/intervals/min_mac${min_mac}"
readonly pheno_dir="data/phenotypes"
readonly in_prefix="ukb_eur_wes_200k"

mkdir -p ${out_dir}

submit_binary_analysis()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/spiros_brava_phenotypes_binary_200k_header.tsv"
  local pheno_file="${pheno_dir}/spiros_brava_phenotypes_binary_200k.tsv"
  local phenotype=$( sed "${SLURM_ARRAY_TASK_ID}q;d" ${pheno_list} )
  submit_intervals "${annotation}" "${phenotype}" "binary" "${pheno_file}"
}

submit_cts_analysis()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local pheno_file="${pheno_dir}/filtered_covar_phenotypes_cts.tsv.gz"
  local phenotype=$( sed "${SLURM_ARRAY_TASK_ID}q;d" ${pheno_list} )
  submit_intervals "${annotation}" "${phenotype}" "cts" "${pheno_file}"
}

submit_intervals()
{
  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local trait=${3?Error: Missing arg3 (trait)}
  local pheno_file=${4?Error: Missing arg4 (pheno_file)}
  local genes="${in_dir}/${in_prefix}_${phenotype}_${annotation}.tsv.gz"
  local out_prefix="${out_dir}/${in_prefix}_${phenotype}_${annotation}"
  if [ ! -z "${phenotype}" ]; then
    if [ -f ${genes} ]; then
      readonly slurm_tasks="${SLURM_ARRAY_TASK_ID}"
      readonly slurm_jname="_filter_genotypes_${1}"
      readonly slurm_lname="_filter_genotypes"
      readonly slurm_project="lindgren.prj"
      readonly slurm_queue="long"
      readonly slurm_shmem="4"
      set -x
      sbatch \
        --account="${slurm_project}" \
        --job-name="${slurm_jname}" \
        --output="${slurm_lname}.log" \
        --error="${slurm_lname}.errors.log" \
        --chdir="${curwd}" \
        --partition="${slurm_queue}" \
        --cpus-per-task="${slurm_shmem}" \
        --array=${slurm_tasks} \
        --parsable \
        "${bash_script}" \
        "${genes}" \
        "${final_sample_list}" \
        "${out_prefix}" \
        "${padding}" \
        "${min_maf}" \
        "${min_info}" \
        "${pheno_file}" \
        "${trait}" \
        "${phenotype}"
    set +x 
  else
      >&2 echo "${genes} (${phenotype}) did not pass significance threshols or does not exist. Skipping.."
    fi
  fi
}

submit_binary_analysis "pLoF_damaging_missense"
#submit_cts_analysis "pLoF_damaging_missense"

#submit_binary_analysis "pLoF"
#submit_cts_analysis "pLoF"

#submit_binary_analysis "damaging_missense"
#submit_cts_analysis "damaging_missense"

#submit_binary_analysis "synonymous"
#submit_cts_analysis "synonymous"






