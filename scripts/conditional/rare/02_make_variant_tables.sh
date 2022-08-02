#!/usr/bin/env bash
#
#$ -N make_variant_tables
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/make_variant_tables.log
#$ -e logs/make_variant_tables.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qa
#$ -t 1-22
#$ -V

set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh

readonly make_script="scripts/conditional/rare/_make_variant_tables.sh"
readonly merge_script="scripts/conditional/rare/_merge_variant_tables.sh"

readonly chr="${SGE_TASK_ID}"
readonly pheno_dir="data/phenotypes"
readonly in_dir="data/conditional/rare/combined"
readonly out_dir="data/conditional/rare/combined/chunks_parallel"

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
readonly queue="short.qc"
readonly nslots=4
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


submit_qc_job() {
  local pheno_list=${1}
  local pheno_file=${2}
  local trait=${3}
  local qsub_main="_pref_c${chr}_${trait}"
  local qsub_mrg="_mrg_c${chr}"
  local pheno_list_csv=$(cat ${pheno_list} | tr "\n" ",")
  local out_prefix_trait="${out_prefix}_${trait}"

  # submit main script
  qsub -N "${qsub_main}" \
      -t ${tasks} \
      -q "${queue}" \
      -pe shmem ${nslots} \
      "${make_script}" \
      "${in_vcf}" \
      "${vcf_lines}" \
      "${lines_per_chunk}" \
      "${chunks}" \
      "${pheno_file}" \
      "${pheno_list_csv}" \
      "${covar_path}" \
      "${out_prefix}"

  # submit merge job
  qsub -N "${qsub_mrg}" \
    -q short.qc@@short.hge \
    -pe shmem 1 \
    -hold_jid "${qsub_main}" \
    "${merge_script}" \
    "${out_prefix}" \
    "${chunks}" \
    "${out_mrg}"

}

submit_binary



