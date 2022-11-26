#!/usr/bin/env bash
#
# @description get perfect LD and MAC by phenotype
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=make_variant_tables
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/make_variant_tables.log
#SBATCH --error=logs/make_variant_tables.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1,17

set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh

readonly curwd=$(pwd)
readonly make_script="scripts/conditional/rare/_make_variant_tables.sh"
readonly merge_script="scripts/conditional/rare/_merge_variant_tables.sh"

readonly chr="${SLURM_ARRAY_TASK_ID}"
readonly pheno_dir="data/phenotypes"
readonly in_dir="data/conditional/rare/combined"
readonly out_dir="data/conditional/rare/combined/chunks"

readonly in_vcf="${in_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense.vcf.bgz"
readonly covar_path="${pheno_dir}/covars1.csv"

readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense"
readonly out_mrg="${out_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense"

mkdir -p ${out_dir}

# integer division with ceiling rounding
int_div () {
  echo $(( $(( ${1} + ${2} - 1)) / ${2} ))
}

# how many chunks should be employed (up to 2 qc cores needed for 1000 lines chunks).
readonly vcf_lines=$( zcat "${in_vcf}" | grep -v "#" | cut -f1 | wc -l )
readonly lines_per_chunk=1000
readonly chunks=$( int_div ${vcf_lines} ${lines_per_chunk})

# hpc paramaters
readonly queue="short"
readonly project="lindgren.prj"
readonly nslots=2
readonly tasks="1-${chunks}"

submit_binary(){
  local trait="binary"
  local pheno_list="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
  local pheno_file="${pheno_dir}/curated_covar_phenotypes_binary_200k.tsv" 
  submit_qc_job ${pheno_list} ${pheno_file} ${trait}
}

submit_cts(){
  local trait="cts"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local pheno_file="${pheno_dir}/curated_covar_phenotypes_cts_200k.tsv" 
  submit_qc_job ${pheno_list} ${pheno_file} ${trait}
}


submit_qc_job() 
{
  # read input
  local pheno_list=${1}
  local pheno_file=${2}
  local trait=${3}
  
  # process input
  local pheno_list_csv=$(cat ${pheno_list} | tr "\n" ",")
  local out_prefix_trait="${out_prefix}_${trait}"

  # need to fix this ....
  local slurm_main_tasks="${tasks}"
  local slurm_main_jname="_pref_c${chr}_${trait}"
  local slurm_main_lname="${out_prefix}"
  local slurm_main_project="${project}"
  local slurm_main_queue="${queue}"
  local slurm_main_nslots="${nslots}"
  readonly make_jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_main_jname}" \
    --output="${slurm_main_lname}.log" \
    --error="${slurm_main_lname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_main_queue}" \
    --cpus-per-task="${slurm_main_nslots}" \
    --array=${slurm_main_tasks} \
    --parsable \
    "${make_script}" \
    "${in_vcf}" \
    "${vcf_lines}" \
    "${lines_per_chunk}" \
    "${chunks}" \
    "${pheno_file}" \
    "${pheno_list_csv}" \
    "${covar_path}" \
    "${out_prefix}" )

  local slurm_merge_jname="_mrg_c${chr}"
  local slurm_merge_lname="${out_prefix}"
  local slurm_merge_project="${project}"
  local slurm_merge_queue="${queue}"
  local slurm_merge_nslots="1"
  readonly merge_jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_merge_jname}" \
    --output="${slurm_merge_lname}.log" \
    --error="${slurm_merge_lname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_merge_queue}" \
    --cpus-per-task="${slurm_merge_nslots}" \
    --dependency="afterok:${make_jid}" \
    --open-mode="append" \
    --parsable \
    "${merge_script}" \
    "${out_prefix}" \
    "${chunks}" \
    "${out_mrg}" )

}

submit_binary



