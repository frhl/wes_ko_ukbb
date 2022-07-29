#!/usr/bin/env bash
#
#$ -N spa_cond_rare
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/spa_cond_rare.log
#$ -e logs/spa_cond_rare.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 2
#$ -tc 10
#$ -V

set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh

readonly rscript="scripts/conditional/rare/_qc_phenotypes.sh"
readonly rcombine="scripts/conditional/rare/_combine_ac.R"

readonly chr="${SGE_TASK_ID}"
readonly pheno_dir="data/phenotypes"
readonly in_dir="data/conditional/rare/combined"
readonly out_dir="data/conditional/rare/combined"

readonly pheno_cts_path="${pheno_dir}/curated_covar_phenotypes_cts_200k.tsv"
readonly pheno_bin_path="${pheno_dir}/curated_covar_phenotypes_binary_200k.tsv"

readonly in_vcf="${in_dir}/ukb_eur_wes_200k_chrCHR_maf0to5e-2_pLoF_damaging_missense.vcf.bgz"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chrCHR_maf0to5e-2_pLoF_damaging_missense_ld"
readonly covar_path="${pheno_dir}/covars1.csv"

readonly phenotypes_cts=$(cat "${pheno_dir}/filtered_phenotypes_cts_manual.tsv" | tr "\n" "," | sed 's/\(.*\),/\1 /' )
readonly phenotypes_bin=$(cat "${pheno_dir}/filtered_phenotypes_binary_header.tsv" | tr "\n" "," | sed 's/\(.*\),/\1 /' )

submit_binary(){
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  local pheno_file=${phenotypes_bin}
  submit_qc_job
}


submit_qc_job() {
  mkdir -p ${step2_dir}
  set -x
  qsub -N "${qsub_spa_name}" \
    -t ${tasks} \
    -q "${queue}" \
    -pe shmem ${nslots} \
    "${rscript}" \
    "${phenotype}" \
    "${pheno_file}" \
    "${covar_path}" \
    "${in_vcf}" \
  set +x
}

readonly tasks=21
readonly queue="short.qc"
submit_binary


