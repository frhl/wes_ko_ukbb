#!/usr/bin/env bash
#
#$ -N prefilter_phenotypes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prefilter_phenotypes.log
#$ -e logs/prefilter_phenotypes.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 21
#$ -tc 10
#$ -V

set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh

readonly rscript="scripts/conditional/rare/_prefilter_phenotypes.sh"

readonly chr="${SGE_TASK_ID}"
readonly pheno_dir="data/phenotypes"
readonly in_dir="data/conditional/rare/combined"
readonly out_dir="data/conditional/rare/combined"

readonly pheno_cts_path="${pheno_dir}/curated_covar_phenotypes_cts_200k.tsv"

readonly in_vcf="${in_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense.vcf.bgz"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense_ld"
readonly covar_path="${pheno_dir}/covars1.csv"

readonly phenotypes_cts="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
readonly phenotypes_bin="${pheno_dir}/filtered_phenotypes_binary_header.tsv"

mkdir -p ${out_dir}

submit_binary(){
  local pheno_list="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
  local pheno_file="${pheno_dir}/curated_covar_phenotypes_binary_200k.tsv"
  submit_qc_job ${pheno_list} ${pheno_file}
}


submit_qc_job() {
  local pheno_list=${1}
  local pheno_file=${2}
  set -x
  qsub -N "_pref_c${chr}" \
    -t ${tasks} \
    -q "${queue}" \
    -pe shmem ${nslots} \
    "${rscript}" \
    "${pheno_list}" \
    "${pheno_file}" \
    "${covar_path}" \
    "${in_vcf}" \
    "${out_prefix}"
  set +x
}

readonly tasks=1-2
readonly queue="short.qc"
readonly nslots=2
submit_binary


