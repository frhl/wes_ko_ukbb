#!/usr/bin/env bash
#
#$ -N fit_simulated_alt
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/fit_simulated_alt.log
#$ -e logs/fit_simulated_alt.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 1-5
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

readonly grm_dir="data/saige/grm/input"
readonly grm_mtx="${grm_dir}/211102_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"

readonly conditioning_markers=""
readonly min_mac=2
readonly chr=22
readonly tasks=${chr}
readonly queue="short.qf"
readonly nslots=1

# wrapper for running main script
run_with_params() {
    submit_spa ${1} ${2} ${3} ${4} ${5} ${6}
}

# main script
submit_spa()
{
  local K=0.1
  local h2=${1}
  local var_beta=${2}
  local var_theta=${3}
  local pi_beta=${4}
  local pi_theta=${5}
  local seed=${6}

  local vars="${var_beta}_${var_theta}"
  local pis="${pi_beta}_${pi_theta}"


  # change "case" to "y" for quant
  local prefix="ukb_eur_h2_${h2}_var_${vars}_pi_${pis}_K${K}_seed${seed}_chr${chr}"
  local saige_prefix="${prefix}_y_${SGE_TASK_ID}"
  
  #local phenotype="y_${SGE_TASK_ID}"
  local phenotype="y_${SGE_TASK_ID}"
  submit_spa_with_csqs "${phenotype}" "${saige_prefix}" "cts"
}


submit_spa_with_csqs()
{
  local phenotype=${1?Error: Missing arg2 (phenotype)}
  local saige_prefix=${2?Error: Missing arg3 (saige_prefix)}
  local trait=${3?Error: Missing arg4 (trait)}
  local step1_dir="data/simulation/saige/step1/${trait}"
  local step2_dir="data/simulation/saige/step2/${trait}"
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
    "${grm_mtx}" \
    "${grm_sam}" \
    "${min_mac}" \
    "${out_prefix}.txt" \
    "${conditioning_markers}"
  set +x
}

readonly annotation="pLoF_damaging_missense"

run_with_params 0.002 0.10 0.10 0.20 0.20 201
run_with_params 0.002 5.00 0.10 0.20 0.20 201
run_with_params 0.002 0.10 5.00 0.20 0.20 201

run_with_params 0.003 0.10 0.10 0.20 0.20 201
run_with_params 0.003 5.00 0.10 0.20 0.20 201
run_with_params 0.003 0.10 5.00 0.20 0.20 201

run_with_params 0.005 0.10 0.10 0.20 0.20 201
run_with_params 0.005 5.00 0.10 0.20 0.20 201
run_with_params 0.005 0.10 5.00 0.20 0.20 201

