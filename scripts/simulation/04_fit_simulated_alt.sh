#!/usr/bin/env bash
#
#$ -N fit_simulated_alt
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/fit_simulated_alt.log
#$ -e logs/fit_simulated_alt.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 1-50
#$ -tc 5
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly vcf_dir="data/simulation/knockouts"
readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"

readonly spa_script="scripts/_spa_test.sh"
readonly merge_script="scripts/_spa_merge.sh"
readonly vcf_prefix="ukb_wes_union_calls_100k_encoded_chr22"

readonly grm_dir="data/saige/grm/input/dnanexus"
readonly grm_mtx="${grm_dir}/ukb_eur_200k_grm_fitted_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"

readonly conditioning_markers=""
readonly min_mac=4
readonly chr=22
readonly tasks=${chr}
readonly queue="short.qe"
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
  local prefix="ukb_wes_union_calls_h2_${h2}_var_${vars}_pi_${pis}_K${K}_seed${seed}_chr${chr}"
  local saige_prefix="${prefix}_case_${SGE_TASK_ID}"
  
  #local phenotype="y_${SGE_TASK_ID}"
  local phenotype="case_${SGE_TASK_ID}"
  submit_spa_with_csqs "${phenotype}" "${saige_prefix}" "binary"
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
  local out_prefix="${step2_dir}/${saige_prefix}"
  local in_vcf="${vcf_dir}/${vcf_prefix}.vcf.bgz"
  submit_spa_job

}


submit_spa_job() {
  >&2 echo "Submitting SPA Job"
  mkdir -p ${step2_dir}
  set -x
  qsub -N "spa_${phenotype}" \
    -o "logs/_fit_simulated_alt.log" \
    -e "logs/_fit_simulated_alt.errors.log" \
    -wd $(pwd) \
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


# this works and results in nice phenotypes
run_with_params 0.00 0.00 0.00 0.01 0.01 100
run_with_params 0.001 0.10 99.0 0.20 0.20 101
run_with_params 0.002 0.10 99.0 0.20 0.20 102
run_with_params 0.005 0.10 99.0 0.20 0.20 103
run_with_params 0.01 0.10 99.0 0.20 0.20 104
run_with_params 0.02 0.10 99.0 0.20 0.20 105
run_with_params 0.05 0.10 99.0 0.20 0.20 106

# simualte additive effects
run_with_params 0.001 99 0.10 0.20 0.20 201
run_with_params 0.002 99 0.10 0.20 0.20 202
run_with_params 0.005 99 0.10 0.20 0.20 203
run_with_params 0.01 99 0.10 0.20 0.20 204
run_with_params 0.02 99 0.10 0.20 0.20 205
run_with_params 0.05 99 0.10 0.20 0.20 206

# this works and results in nice phenotypes
#run_with_params 0.00 0.00 0.00 0.01 0.01 100
#run_with_params 0.001 0.10 99.0 0.20 0.20 101
#run_with_params 0.002 0.10 99.0 0.20 0.20 102
#run_with_params 0.005 0.10 99.0 0.20 0.20 103
#run_with_params 0.01 0.10 99.0 0.20 0.20 104
#run_with_params 0.02 0.10 99.0 0.20 0.20 105
#run_with_params 0.05 0.10 99.0 0.20 0.20 106


run_with_params 0.001 99.0 0.10 0.20 0.20 101
run_with_params 0.002 99.0 0.10 0.20 0.20 102
run_with_params 0.005 99.0 0.10 0.20 0.20 103
run_with_params 0.01 99.0 0.10 0.20 0.20 104
run_with_params 0.02 99.0 0.10 0.20 0.20 105
run_with_params 0.05 99.0 0.10 0.20 0.20 106


#run_with_params 0.001 0.10 99.0 1.00 1.00 601
#run_with_params 0.002 0.10 99.0 1.00 1.00 602
#run_with_params 0.005 0.10 99.0 1.00 1.00 603
#run_with_params 0.01 0.10 99.0 1.00 1.00 604
#run_with_params 0.02 0.10 99.0 1.00 1.00 605
#run_with_params 0.05 0.10 99.0 1.00 1.00 606


#run_with_params 0.001 10.0 0.10 0.20 0.20 601
#run_with_params 0.001 0.10 0.10 0.20 0.20 601
#run_with_params 0.001 0.10 1.00 0.20 0.20 602
#run_with_params 0.001 0.10 10.0 0.20 0.20 603
#run_with_params 0.001 0.10 99.0 0.20 0.20 604

#run_with_params 0.005 10.0 0.10 0.20 0.20 501
#run_with_params 0.005 0.10 0.10 0.20 0.20 501
#run_with_params 0.005 0.10 1.00 0.20 0.20 502
#run_with_params 0.005 0.10 10.0 0.20 0.20 503
#run_with_params 0.005 0.10 99.0 0.20 0.20 504

#run_with_params 0.01 10.0 0.10 0.20 0.20 501
#run_with_params 0.01 0.10 0.10 0.20 0.20 501
#run_with_params 0.01 0.10 1.00 0.20 0.20 502
#run_with_params 0.01 0.10 10.0 0.20 0.20 503
#run_with_params 0.01 0.10 99.0 0.20 0.20 504

#run_with_params 0.10 10.0 0.10 0.20 0.20 501
#run_with_params 0.10 0.10 0.10 0.20 0.20 501
#run_with_params 0.10 0.10 1.00 0.20 0.20 502
#run_with_params 0.10 0.10 10.0 0.20 0.20 503
#run_with_params 0.10 0.10 99.0 0.20 0.20 504

