#!/usr/bin/env bash
#
#$ -N array_spa
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/array_spa.log
#$ -e logs/array_spa.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 1
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"
readonly bash_script="scripts/permute/_submit_chr_spa.sh"

readonly in_dir="data/permute/genes/chrCHR"
readonly out_dir="data/permute/spa/chrCHR"

readonly in_prefix="${in_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_permuted_chrCHR_GENE"
readonly out_prefix="${out_dir}/test_chrCHR_GENE"

readonly overview="data/permute/overview/overview.tsv.gz"

readonly min_mac=5
readonly tasks=21
readonly queue="short.qf"
readonly nslots=1

submit_spa_set_binary()
{
  local pheno_list="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_spa_pair "${annotation}" "${phenotype}" "binary"
}

submit_spa_set_cts()
{
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_spa_pair  "${phenotype}" "cts"
}

submit_spa_pair()
{

  local phenotype=${1?Error: Missing arg2 (phenotype)}
  local trait=${2?Error: Missing arg3 (trait)}

  local step1_dir="data/saige/output/${trait}/step1"
  local step2_dir="data/saige/output/set/${trait}/step2"
  local in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
  local in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
  local out_prefix_pheno="${in_prefix}_${phenotype}"
  submit_spa_gene_job
}


submit_spa_gene_job()
{
  mkdir -p ${step2_dir}
  set -x
  qsub -N "spa_${phenotype}" \
    -o "logs/submit_chr_spa.log" \
    -e "logs/submit_chr_spa.errors.log" \
    -t ${tasks} \
    -q "${queue}" \
    -pe shmem ${nslots} \
    "${bash_script}" \
    "${phenotype}" \
    "${in_prefix}" \
    "${in_gmat}" \
    "${in_var}" \
    "${min_mac}" \
    "${overview}" \
    "${out_prefix_pheno}"
  set +x
}

submit_spa_set_cts

