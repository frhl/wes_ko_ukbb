#!/usr/bin/env bash
#
#$ -N fit_simulated_alt
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/fit_simulated_alt.log
#$ -e logs/fit_simulated_alt.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 1-25
#$ -tc 5
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
readonly min_mac=2
readonly chr=21
readonly tasks=${chr}
readonly queue="short.qf"
readonly nslots=1

submit_spa()
{
  local K=0.1
  local h2_nc=${1}
  local h2_co=${2}
  local h2_ko=${3}
  local pi_nc=${4}
  local pi_co=${5}
  local pi_ko=${6}
  local alpha=${7}
  local beta=${8}
  local theta=${9}

  local h2s="${h2_nc}_${h2_co}_${h2_ko}"
  local pis="${pi_nc}_${pi_co}_${pi_ko}"
  local effects="a${alpha}_b${beta}_t${theta}"

  local prefix="ukb_eur_h2_${h2s}_pi_${pis}_K${K}_${effects}_chr${chr}"
  local saige_prefix="${prefix}_y_cts_${SGE_TASK_ID}"
  local phenotype="y_cts_${SGE_TASK_ID}"
  submit_spa_with_csqs "${phenotype}" "${saige_prefix}" "cts"
}


submit_spa_with_csqs()
{
  local phenotype=${1?Error: Missing arg2 (phenotype)}
  local saige_prefix=${2?Error: Missing arg3 (saige_prefix)}
  local trait=${3?Error: Missing arg4 (trait)}
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

readonly annotation="pLoF_damaging_missense"

# Absence of CH effects
submit_spa 0.00 0.00 0.00 0.00 0.00 0.00 NA NA NA
submit_spa 0.10 0.00 0.00 0.00 0.00 0.00 NA NA NA
submit_spa 0.00 0.10 0.00 0.00 0.00 0.00 NA NA NA
submit_spa 0.10 0.10 0.00 0.00 0.00 0.00 NA NA NA

# simulate CH effects
submit_spa 0.00 0.00 0.02 0.00 0.00 0.10 NA NA NA
submit_spa 0.00 0.10 0.02 0.00 0.00 0.10 NA NA NA
submit_spa 0.00 0.10 0.02 0.00 0.00 0.10 NA NA NA

# simulate effects with thetas
submit_spa 0.00 0.00 0.00 0.00 0.00 0.10 NA NA 0.01
submit_spa 0.00 0.00 0.00 0.00 0.00 0.10 NA NA 0.10


