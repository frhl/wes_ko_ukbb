#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=fit_simulated_alt
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/fit_simulated_alt.log
#SBATCH --error=logs/fit_simulated_alt.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1
#SBATCH --array=1-10

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

# array idx indicate phenotype idx
readonly array_idx=$( get_array_task_id )
readonly chr=22

readonly vcf_dir="data/simulation/knockouts"
readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"

readonly spa_script="scripts/_spa_test.sh"
#readonly merge_script="scripts/_spa_merge.sh"
readonly vcf_prefix="ukb_wes_union_calls_100k_encoded_chr${chr}"

readonly grm_dir="data/saige/grm/input/dnanexus"
readonly grm_mtx="${grm_dir}/ukb_eur_200k_grm_fitted_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"

readonly conditioning_markers=""

# wrapper for running main script
run_with_params() {
    submit_spa ${1} ${2} ${3} ${4} ${5}
}

# main script
submit_spa()
{
  local K=${1}
  local h2=${2}
  local b=${3}
  local pi=${4}
  local seed=${5}


  # change "case" to "y" for quant
  local prefix="ukb_wes_union_calls_h2_${h2}_b_${b}_pi_${pi}_K${K}_seed${seed}_chr${chr}"
  local saige_prefix="${prefix}_case_${array_idx}"

  #local phenotype="y_${array_idx}"
  local phenotype="case_${array_idx}"
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
  local jname="_fit_simulated_alt"
  local lname="logs/_fit_simulated_alt"
  set -x
  local sim_job=$( sbatch \
    --account="${project}" \
    --job-name="${jname}" \
    --output="${lname}.log" \
    --error="${lname}.errors.log" \
    --chdir="$(pwd)" \
    --partition="${queue}" \
    --cpus-per-task="${nslots}" \
    --array="${chr}" \
    --parsable \
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
    "${conditioning_markers}")
}

readonly min_mac=4
readonly project="lindgren.prj"
readonly queue="short"
readonly nslots=1

run_with_params 0.2 0.00 0.5 0.20 101
run_with_params 0.2 0.01 0.5 0.20 101
run_with_params 0.2 0.02 0.5 0.20 101
run_with_params 0.2 0.05 0.5 0.20 101

run_with_params 0.2 0.01 0.0 0.20 101
run_with_params 0.2 0.01 0.5 0.20 101
run_with_params 0.2 0.01 1.0 0.20 101
run_with_params 0.2 0.01 10.0 0.20 101

