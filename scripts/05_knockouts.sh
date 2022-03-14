#!/usr/bin/env bash
#
#
#$ -N knockouts
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/test.log
#$ -e logs/test.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly in_dir="data/mt/csqs"
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/knockouts/alt_filtered"

readonly knockout_script="scripts/_knockouts.sh"
readonly exclude="data/genes/220310_common_plofs_to_exclude.txt"
readonly input_path="${in_dir}/ukb_eur_wes_200k_annot_chrCHR.mt"
readonly input_type="mt"

readonly out_prefix="${out_dir}/ukb_eur_wes_200k"
readonly out_type="vcf"

# Note: ~20 slots are needed for running chr1 
# when using aggr_method="collect" on short.qe
readonly tasks="1-22"
readonly queue="short.qe"
readonly nslots=20

readonly phase=""
readonly seed=""
readonly only_vcf=""
readonly aggr_method="collect" # either fasts or collect

readonly maf_lb=0
readonly maf_ub=5e-2

mkdir -p ${out_dir}

submit_knockout_job() 
{
  local prefix="${out_prefix}_chrCHR_maf${maf_lb}to${maf_ub}${sex:+_${sex}}"
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
    "${prefix}" \
    "${out_type}" \
    "${aggr_method}" \
    "${phase}" \
    "${seed}" \
    "${only_vcf}" \
  set +x
}

#submit_knockout_job "0" "5e-2" "" "pLoF"
#submit_knockout_job "0" "5e-2" "" "damaging_missense"
submit_knockout_job "0" "5e-2" "" "pLoF_damaging_missense"
#submit_knockout_job "0" "5e-2" "" "pLoF,LC,damaging_missense"

#submit_knockout_job 0 0.05 "" "pLoF"
#submit_knockout_job 0 0.05 "" "synonymous"
#§submit_knockout_job 0 0.05 "" "ptv,LC"
