#!/usr/bin/env bash
#
#$ -N array_spa
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/array_spa.log
#$ -e logs/array_spa.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 44
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"
readonly bash_script="scripts/permute/_array_spa.sh"

readonly in_dir="data/permute/permutations/chrCHR"
readonly out_dir="data/permute/spa/chrCHR"
readonly overview_dir="data/permute/overview"

readonly in_prefix="${in_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_permuted_chrCHR_GENE"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_permuted_chrCHR_GENE"

readonly overview="${overview_dir}/overview.tsv.gz"
readonly gene_spa="${overview_dir}/overview_genes.tsv.gz"

readonly min_mac=5
readonly tasks=22
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
  local out_prefix_pheno="${out_prefix}_${phenotype}"
  submit_spa_gene_job
}


submit_spa_gene_job()
{
  mkdir -p ${step2_dir}
  set -x
  qsub -N "spa_${phenotype}" \
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
    "${gene_spa}" \
    "${out_prefix_pheno}" \
    "${out_dir}"
  set +x
}

submit_spa_set_cts

