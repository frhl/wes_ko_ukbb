#!/usr/bin/env bash
#
#
#$ -N rand_phase_ko
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/rand_phase_ko.log
#$ -e logs/rand_phase_ko.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly in_dir="data/mt"
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/knockouts/null"

readonly knockout_script="scripts/_phased_knockouts.sh"
readonly input_path="${in_dir}/ukb_wes_200k_merged_chrCHR.mt"
readonly input_type="mt"

readonly out_type="vcf"

readonly af_min=""
readonly af_max=""
readonly phase="random"
readonly only_vcf="yes"

readonly tasks="22"
readonly queue="short.qe"
readonly nslots=2

mkdir -p ${out_dir}

submit_knockout_random_job() 
{
  local maf_lb=${1}
  local maf_ub=${2}
  local sex=${3}
  local csq=${4}
  local prefix="${out_prefix}_chrCHR_"maf${maf_lb}_${maf_ub}${sex:+_${sex}}
  local qsub_name=$( echo ${csq} | tr "," "_")
  
  set -x
  qsub -N "_${qsub_name}" \
    -o "logs/_knockouts_random.log" \
    -e "logs/_knockouts_random.errors.log" \
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
    "${phase}" \
    "${only_vcf}"
  set +x
}


#for i in $(seq 1 4); do 
  out_prefix="${out_dir}/test_varid_ukb_wes_200k_rand${i}"
  submit_knockout_random_job "0" "5e-2" "" "pLoF,damaging_missense"
#done

#submit_knockout_random_job 0 0.05 "" "pLoF"
#submit_knockout_random_job 0 0.05 "" "synonymous"
#submit_knockout_random_job 0 0.05 "" "ptv,LC"
