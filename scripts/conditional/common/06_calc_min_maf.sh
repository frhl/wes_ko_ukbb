#!/usr/bin/env bash
#
# extract genotypes in regions near genes that are significant in primary analysis
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=calc_min_maf
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/calc_min_maf.log
#SBATCH --error=logs/calc_min_maf.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-350


set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )

readonly bash_script="scripts/conditional/common/_calc_min_maf.sh"
readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'

readonly min_mac=4
readonly in_dir="data/conditional/common/intervals/2024/min_mac${min_mac}"
readonly out_dir="data/conditional/common/intervals/2024/min_mac${min_mac}"
readonly pheno_dir="data/phenotypes"
readonly in_prefix="ukb_eur_wes_200k"

mkdir -p ${out_dir}

submit_binary_analysis()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/dec22_phenotypes_binary_200k_header.tsv"
  local pheno_file="${pheno_dir}/dec22_phenotypes_binary_200k.tsv.gz"
  local phenotype=$( sed "${task_id}q;d" ${pheno_list} )
  submit_intervals "${annotation}" "${phenotype}" "binary" "${pheno_file}"
}

submit_cts_analysis()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local pheno_file="${pheno_dir}/filtered_covar_phenotypes_cts.tsv.gz"
  local phenotype=$( sed "${task_id}q;d" ${pheno_list} )
  submit_intervals "${annotation}" "${phenotype}" "cts" "${pheno_file}"
}

submit_intervals()
{
  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local trait=${3?Error: Missing arg3 (trait)}
  local pheno_file=${4?Error: Missing arg4 (pheno_file)}
  local intervals="${in_dir}/${in_prefix}_${phenotype}_${annotation}.ht"
  local out_prefix="${out_dir}/${in_prefix}_${phenotype}_${annotation}"
  if [ ! -z "${phenotype}" ]; then
    if [ -d "${intervals}" ]; then
      readonly slurm_jname="_calc_min_maf_${1}"
      readonly slurm_lname="logs/_calc_min_maf"
      readonly slurm_project="lindgren.prj"
      readonly slurm_queue="short"
      readonly sge_queue="short.qc"
      readonly slurm_shmem="1"
      if [ "${cluster}" = "slurm" ]; then
        sbatch \
          --account="${slurm_project}" \
          --job-name="${slurm_jname}" \
          --output="${slurm_lname}.log" \
          --error="${slurm_lname}.errors.log" \
          --chdir="$(pwd)" \
          --partition="${slurm_queue}" \
          --cpus-per-task="${slurm_shmem}" \
          --array=${task_id} \
          --parsable \
          "${bash_script}" \
          "${intervals}" \
          "${out_prefix}" \
          "${pheno_file}" \
          "${trait}" \
          "${phenotype}"
      elif [ "${cluster}" = "sge" ]; then
        qsub -N "${slurm_jname}" \
          -o "${slurm_lname}.log" \
          -e "${slurm_lname}.errors.log" \
          -P lindgren.prjc \
          -wd $(pwd) \
          -t ${task_id} \
          -q "${sge_queue}" \
          -pe shmem ${slurm_shmem} \
          "${bash_script}" \
          "${intervals}" \
          "${out_prefix}" \
          "${pheno_file}" \
          "${trait}" \
          "${phenotype}"
      fi 
    else
      >&2 echo "${intervals} Does not exist. Skipping.."
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






