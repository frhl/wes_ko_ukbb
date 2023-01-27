#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=co_occurence
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/co_occurence.log
#SBATCH --error=logs/co_occurence.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-2
#
#$ -N co_occurence
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/co_occurence.log
#$ -e logs/co_occurence.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 51-300
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly curwd=$(pwd)
readonly bash_script="scripts/_co_occurence.sh"
readonly rscript="scripts/_co_occurence.R"

readonly pheno_dir="data/phenotypes"
readonly phenotypes="${pheno_dir}/dec22_phenotypes_binary_200k.tsv.gz"
readonly pheno_list="${pheno_dir}/dec22_phenotypes_binary_200k_header.tsv"

readonly out_dir="data/test/"
readonly out_prefix="${out_dir}/ukb_wes_200k"
readonly annotation="pLoF_damaging_missense"

readonly cluster=$( get_current_cluster )
readonly index=$( get_array_task_id )
readonly phenotype=$( sed "${index}q;d" ${pheno_list} )

readonly slurm_tasks="1-22" 
readonly slurm_jname="_co_occurence_${phenotype}"
readonly slurm_lname="logs/_co_occurence"
readonly slurm_project="lindgren.prj"
readonly slurm_queue="short"
readonly sge_queue="short.qc"
readonly slurm_nslots="1"

if [ ${cluster} == "slurm" ]; then
  sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname}.log" \
    --error="${slurm_lname}.errors.log" \
    --chdir="$(pwd)" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --array=${slurm_tasks} \
    --parsable \
    "${bash_script}" \
    "${pheno_file}" \
    "${phenotype}" \
    "${annotation}" \
    "${out_prefix}"
elif [ "${cluster}" == "sge" ]; then
  qsub -N "${slurm_jname}" \
    -q ${sge_queue} \
    -pe shmem ${slurm_nslots} \
    -P lindgren.prjc \
    -o "${slurm_lname}.log" \
    -e "${slurm_lname}.errors.log" \
    -t "${slurm_tasks}" \
    -wd $(pwd) \
    "${bash_script}" \
    "${pheno_file}" \
    "${phenotype}"
    "${annotation}" \
    "${out_prefix}"
fi





