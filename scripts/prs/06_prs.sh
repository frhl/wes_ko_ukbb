#!/usr/bin/env bash
#
#$ -N prs
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prs.log
#$ -e logs/prs.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 1-2
#$ -tc 1
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly bash_script="scripts/prs/_prs.sh"
readonly rscript="scripts/prs/06_prs.R"

readonly ldsc_dir="data/prs/ldsc"
readonly pred_dir="data/prs/hapmap/ukb_500k"
readonly ld_dir="data/prs/hapmap/ld/matrix"
readonly pheno_dir="data/phenotypes"
readonly out_dir="data/prs/scores/auto"

readonly index=${SGE_TASK_ID}

readonly file_cts="${pheno_dir}/filtered_phenotypes_cts.tsv"
readonly pheno_list_cts="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
readonly phenotype_cts=$( sed "${index}q;d" ${pheno_list_cts} )

readonly file_binary="${pheno_dir}/filtered_phenotypes_binary.tsv"
readonly pheno_list_binary="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
readonly phenotype_binary=$( sed "${index}q;d" ${pheno_list_binary} )

mkdir -p ${out_dir}

submit_ldpred2()
{ 
  local method=${1}
  local cores=${2}
  local phenotype=${3}
  local pred="${pred_dir}/ukb_hapmap_500k_eur_chrCHR.bed"
  local ldsc="${ldsc_dir}/ldsc_${phenotype}.rds"
  local out_prefix="${out_dir}/prs_${method}_${phenotype}_chrCHR"
  set -x
  qsub -N "_prs_${phenotype}" \
    -t 1-22 \
    -q short.qc@@short.hga \
    -pe shmem "${cores}" \
    "${bash_script}" \
    "${rscript}" \
    "${pred}" \
    "${ldsc}" \
    "${ld_dir}" \
    "${method}" \
    "${out_prefix}"
  set +x 
}


submit_ldpred2 "auto" "6" "${phenotype_cts}"
submit_ldpred2 "auto" "6" "${phenotype_binary}"

