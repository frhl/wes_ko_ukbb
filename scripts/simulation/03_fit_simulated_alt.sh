#!/usr/bin/env bash
#
#$ -N fit_simulated_alt
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/fit_simulated_alt.log
#$ -e logs/fit_simulated_alt.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 1
#$ -tc 1
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly vcf_dir="data/knockouts/alt"
readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"

readonly spa_script="scripts/_spa_test.sh"
readonly merge_script="scripts/_spa_merge.sh"
readonly vcf_prefix="ukb_eur_wes_200k_chrCHR_maf0to5e-2"

readonly conditioning_markers=""
readonly min_mac=1
readonly chr=21
readonly tasks=${chr}
readonly queue="short.qf"
readonly nslots=1

submit_spa_cts_with_csqs()
{
  local saige_prefix="${1?Error: Missing arg1 (saige_prefix)}"
  local annotation="${2?Error: Missing arg2 (annotation)}"
  local phenotype="y_cts_${SGE_TASK_ID}"
  submit_spa_with_csqs "${annotation}" "${phenotype}" "${saige_prefix}" "cts"
}


submit_spa_with_csqs()
{
  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local saige_prefix=${3?Error: Missing arg3 (saige_prefix)}
  local trait=${4?Error: Missing arg4 (trait)}
  local step1_dir="data/simulation/saige/step1"
  local step2_dir="data/simulation/saige/step2"
  local in_gmat="${step1_dir}/${saige_prefix}.rda"
  local in_var="${step1_dir}/${saige_prefix}.varianceRatio.txt"
  local out_prefix="${step2_dir}/${saige_prefix}_${annotation}"
  local in_vcf="${vcf_dir}/${vcf_prefix}_${annotation}.vcf.bgz"
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
    "${out_prefix}.txt" \
    "${conditioning_markers}"
  set +x
}

submit_spa_cts_with_csqs "ukb_eur_h2_0_0_pi_NA_NA_K_1e-1_chr21_y_cts_${SGE_TASK_ID}" "pLoF_damaging_missense"
#submit_spa_cts_with_csqs "ukb_eur_h2_0_pi_0_K_1e-1_chr21_cts${SGE_TASK_ID}" "pLoF_damaging_missense"



