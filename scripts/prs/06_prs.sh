#!/usr/bin/env bash
#
#$ -N prs
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prs.log
#$ -e logs/prs.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 1-75
#$ -tc 1
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly bash_script="scripts/prs/_prs.sh"
readonly rscript="scripts/prs/06_prs.R"
readonly clean_script="scripts/prs/_prs_clean.sh"
readonly aggr_script="scripts/prs/_prs_aggr.sh"

readonly ldsc_dir="data/prs/ldsc"
readonly pred_dir="data/prs/hapmap/ukb_500k/validation"
readonly ld_dir="data/prs/hapmap/ld/matrix_unrel_kin"
readonly pheno_dir="data/phenotypes"
readonly out_dir="data/prs/scores/auto"
readonly mrg_dir="data/prs/scores"

readonly index=${SGE_TASK_ID}

readonly file_cts="${pheno_dir}/curated_covar_phenotypes_cts.tsv.gz"
readonly pheno_list_cts="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
readonly phenotype_cts=$( sed "${index}q;d" ${pheno_list_cts} )

readonly file_binary="${pheno_dir}/curated_covar_phenotypes_binary.tsv.gz"
readonly pheno_list_binary="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
readonly phenotype_binary=$( sed "${index}q;d" ${pheno_list_binary} )
readonly impute="mean2"

mkdir -p ${out_dir}

submit_ldpred2()
{ 
  local method=${1}
  local cores=${2}
  local phenotype=${3}
  local pred="${pred_dir}/ukb_hapmap_500k_eur_chrCHR.bed"
  local ldsc="${ldsc_dir}/ldsc_${phenotype}.rds"
  local out_prefix="${out_dir}/prs_${method}_${phenotype}_chrCHR"

  local qsub_fit="_prs_${phenotype}"
  local qsub_aggr="_aggr_${phenotype}"
  local qsub_clean="_clean_${phenotype}"

  if [ ! -z ${phenotype} ]; then
    # fit actual pgs
    fit_pgs
    # aggregate into matrix
    aggr_pgs
    # remove disk backing files
    clean_pgs
  fi

}


fit_pgs()
{
  set -x
  qsub -N "${qsub_fit}" \
    -t ${tasks} \
    -q short.qc@@short.hga \
    -pe shmem "${cores}" \
    "${bash_script}" \
    "${rscript}" \
    "${pred}" \
    "${ldsc}" \
    "${ld_dir}" \
    "${method}" \
    "${impute}" \
    "${out_prefix}"
  set +x
}


aggr_pgs()
{
  set -x
  qsub -N "${qsub_aggr}" \
    -q test.qc \
    -pe shmem 1 \
    -hold_jid "_prs_${phenotype}" \
    "${aggr_script}" \
    "${phenotype}" \
    "${out_dir}" \
    "${mrg_dir}"
  set +x

}


clean_pgs()
{
  set -x
  qsub -N "${qsub_clean}" \
    -t ${tasks} \
    -q test.qc \
    -pe shmem 1 \
    -hold_jid_ad "_prs_${phenotype}" \
    "${clean_script}" \
    "${pred}" \
    "${out_prefix}"
  set +x 
}

readonly tasks=1-22
#submit_ldpred2 "auto" "6" "${phenotype_cts}_int"
submit_ldpred2 "auto" "6" "${phenotype_binary}"

