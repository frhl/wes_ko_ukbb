#!/usr/bin/env bash
#
#$ -N absence_of_effect_null
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/absence_of_effect_null.log
#$ -e logs/absence_of_effect_null.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc@@short.hge
#$ -t 0
#$ -V

set -o errexit
set -o nounset

set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh

readonly vcf_dir="data/knockouts/alt"
readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"

readonly spa_script="scripts/_spa_test.sh"
readonly merge_script="scripts/_spa_merge.sh"
readonly in_prefix="ukb_eur_wes_200k_chrCHR"

readonly conditioning_markers=""
readonly min_mac=4
readonly tasks=1-22
readonly queue="short.qf"
readonly nslots=1

submit_spa_cts_with_csqs()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype="y${SGE_TASK_ID}"
  submit_spa_with_csqs "${annotation}" "${phenotype}" "cts"
}


submit_spa_with_csqs()
{
  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local trait=${3?Error: Missing arg3 (trait)}
  local step1_dir="data/saige/output/${trait}/step1"
  local step2_dir="data/saige/output/${trait}/step2"
  local in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
  local in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
  local out_prefix="${step2_dir}/${in_prefix}_${maf}_${phenotype}_${annotation}"
  local in_vcf="${vcf_dir}/${in_prefix}_${maf}_${annotation}.vcf.bgz"
  submit_spa_job

}


submit_spa_job() {
  mkdir -p ${step2_dir}
  set -x
  qsub -N "spa_${phenotype}_${annotation}" \
    -t ${tasks} \
    -q "${queue}" \
    -pe shmem ${nslots} \
    "${spa_script}" \
    "${phenotype}" \
    "${in_vcf}" \
    "${in_vcf}.csi" \
    "${in_gmat}" \
    "${in_var}" \
    "${min_mac}" \
    "${out_prefix}" \
    "${conditioning_markers}"
  set +x
}


