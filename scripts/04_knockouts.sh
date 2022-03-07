#!/usr/bin/env bash
#
#
#$ -N knockouts
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/knockouts.log
#$ -e logs/knockouts.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly in_dir="data/mt/annotated"
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/knockouts/alt/new"

readonly knockout_script="scripts/_knockouts.sh"
readonly input_path="${in_dir}/ukb_eur_wes_200k_annot_chrCHR.mt"
readonly input_type="mt"

readonly af_min=""
readonly af_max=""
readonly phase=""

readonly out_prefix="${out_dir}/ukb_eur_wes_200k"
readonly out_type="vcf"

readonly tasks="1-22"
readonly queue="short.qe"
readonly nslots=4

mkdir -p ${out_dir}

submit_knockout_job() 
{
  local maf_lb=${1}
  local maf_ub=${2}
  local sex=${3}
  local csq=${4}
  local prefix="${out_prefix}_chrCHR_maf${maf_lb}to${maf_ub}_${sex:+_${sex}}"
  local qsub_name=$( echo ${csq} | tr "," "_")
  
  set -x
  qsub -N "_${qsub_name}" \
    -o "logs/_knockouts.log" \
    -e "logs/_knockouts.errors.log" \
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
    "${out_type}" \
    "${phase}"
  set +x
}

submit_knockout_job "0" "5e-2" "" "pLoF"
submit_knockout_job "0" "5e-2" "" "damaging_missense"
submit_knockout_job "0" "5e-2" "" "pLoF,damaging_missense"
submit_knockout_job "0" "5e-2" "" "pLoF,LC,damaging_missense"

#submit_knockout_job 0 0.05 "" "pLoF"
#submit_knockout_job 0 0.05 "" "synonymous"
#submit_knockout_job 0 0.05 "" "ptv,LC"
