#!/usr/bin/env bash
#
#
#$ -N knockouts_rand_phase
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/knockouts_rand_phase.log
#$ -e logs/knockouts_rand_phase.errors.log
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

readonly knockout_script="scripts/_knockouts.sh"
readonly exclude="data/genes/220310_common_plofs_to_exclude.txt"
readonly input_path="${in_dir}/ukb_eur_wes_200k_annot_chrCHR.mt"
readonly input_type="mt"

readonly af_min=""
readonly af_max=""

readonly only_vcf=""
readonly checkpoint=""
readonly aggr_method="fast"
readonly phase="random"

readonly tasks="1-22"
readonly nslots=4
readonly queue="short.qa"

submit_knockout_random_job() 
{
  local maf_lb=${1}
  local maf_ub=${2}
  local sex=${3}
  local csq=${4}
  local seed=${5}
  local prefix="${out_prefix}_chrCHR_maf${maf_lb}to${maf_ub}${sex:+_${sex}}"
  local qsub_name=$( echo ${csq} | tr "," "_")
  
  set -x
  qsub -N "_${qsub_name}" \
    -o "logs/_knockouts_rand_phase.log" \
    -e "logs/_knockouts_rand_phase.errors.log" \
    -t ${tasks} \
    -q "${queue}" \
    -pe shmem ${nslots} \
    "${knockout_script}" \
    "${input_path}" \
    "${input_type}" \
    "${maf_lb}" \
    "${maf_ub}" \
    "${sex}" \
    "${csq}" \
    "${prefix}" \
    "${out_type}" \
    "${aggr_method}" \
    "${phase}" \
    "${seed}" \
    "${only_vcf}" \
    "${exclude}" 
  set +x
}


out_type="vcf" 
for seed in $(seq 6 6); do
  out_dir="data/knockouts/null/seed${seed}"
  out_prefix="${out_dir}/ukb_eur_wes_200k_rand_phase"
  mkdir -p ${out_dir}
  submit_knockout_random_job "0" "5e-2" "" "pLoF,damaging_missense" "${seed}"
done



