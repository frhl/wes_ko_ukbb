#!/usr/bin/env bash
#
# probabilistic model for human knockouts
#
#$ -N knockouts
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/knockouts.log
#$ -e logs/knockouts.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qf
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly in_dir="data/mt"
readonly spark_dir="data/tmp/spark"
readonly out_dir="derived/knockouts/220207"

readonly knockout_script="scripts/_knockouts.sh"
readonly input_path="${in_dir}/ukb_wes_200k_merged_chrCHR.mt"
readonly input_type="mt"

readonly af_min=""
readonly af_max=""

readonly out_prefix="${out_dir}/test_varid_ukb_wes_200k"
readonly out_type="vcf"

readonly tasks="21"
readonly queue="short.qe"
readonly nslots=3


submit_knockout_job() 
{
  mkdir -p ${out_dir}
  local maf_lb=${1}
  local maf_ub=${2}
  local sex=${3}
  local csq=${4}
  local prefix="${out_prefix}_maf${maf_lb}_${maf_ub}_${sex:+_${sex}}chrCHR"
  local qsub_name=$( echo ${csq} | tr "," "_")
  
  set -x
  qsub -N "_${qsub_name}" \
    -t ${tasks} \
    -q "${queue}" \
    -pe shmem ${nslots} \
    "${knockout_script}" \
    "${input_path}" \
    "${input_type}" \
    "${af_min}" \
    "${af_max}" \
    "${maf_lb}" \
    "${maf_ub}" \
    "${sex}" \
    "${csq}" \
    "${prefix}" \
    "${out_type}"
  set +x
}

submit_knockout_job "0" "1e-2" "" "pLoF,damaging_missense"
#submit_knockout_job 0 0.01 "" "pLoF"
#submit_knockout_job 0 0.01 "" "synonymous"
#submit_knockout_job 0 0.01 "" "ptv,LC"
