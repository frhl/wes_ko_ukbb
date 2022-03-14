#!/usr/bin/env bash
#
#$ -N filter_csqs
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/filter_csqs.log
#$ -e logs/filter_csqs.errors.log
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
readonly out_dir="data/mt/csqs"

readonly bash_script="scripts/_filter_csqs.sh"
readonly exclude="data/genes/220310_common_plofs_to_exclude.txt"
readonly input_path="${in_dir}/ukb_eur_wes_200k_annot_chrCHR.mt"
readonly input_type="mt"

readonly out_prefix="${out_dir}/ukb_eur_wes_200k"
readonly out_type="vcf"

readonly af_min=""
readonly af_max=""

readonly queue="short.qe"
readonly tasks=21
readonly nslots=2

mkdir -p ${out_dir}

submit_filter_job() 
{
  local maf_lb=${1}
  local maf_ub=${2}
  local sex=${3}
  local csq=${4}
  local csq_name="$( echo $csq | tr ',' '_' )"
  local prefix="${out_prefix}_chrCHR_maf${maf_lb}to${maf_ub}${sex:+_${sex}}_${csq_name}"
  
  set -x
  qsub -N "_filter_${csq}" \
    -o "logs/_filter_csqs.log" \
    -e "logs/_filter_csqs.errors.log" \
    -t "${tasks}" \
    -q "${queue}" \
    -pe shmem ${nslots} \
    "${bash_script}" \
    "${input_path}" \
    "${input_type}" \
    "${maf_lb}" \
    "${maf_ub}" \
    "${sex}" \
    "${csq}" \
    "${prefix}" \
    "${out_type}" \
    "${exclude}"
  set +x
}

# note, we need to use comma here, as we can't split on '_'
submit_filter_job "0" "5e-2" "" "pLoF,damaging_missense"

#submit_knockout_job 0 0.05 "" "synonymous"


