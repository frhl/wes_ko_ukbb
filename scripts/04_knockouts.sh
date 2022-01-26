#!/usr/bin/env bash
#
# Combine phased and unphased data to make a probabilistic model
# for human knockouts.
#
# 1) A fake VCF with probablistic encoding of KO (for SAIGE+ input)
# 2) A matrix containg the variant consequences in each gene for each sample
# 3) A ko probabiltiy matrix containg individuals, genes and probability of KO
# 4) A ko probability matrix with the above + variants involved in KO.
#
#$ -N knockout
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
readonly out_dir="derived/knockouts/220126"

readonly knockout_script="scripts/_knockouts.sh"
readonly input_path="${in_dir}/ukb_wes_200k_merged_chrCHR.mt"
readonly input_type="mt"

readonly af_min=""
readonly af_max=""

readonly out_prefix="${out_dir}/ukb_wes_200k"
readonly out_type="vcf"

submit_knockout_job() 
{
  
  local maf_lb=${1}
  local maf_ub=${2}
  local sex=${3}
  local csq=${4}

  local prefix="${out_prefix}_maf${maf_lb}_${maf_ub}_${sex:+_${sex}}chrCHR"
  local qsub_name=$( echo ${csq} | tr "," "_")
  
  set -x
  qsub -N "_ko_${qsub_name}" \
    -t 21 \
    -q "short.qa" \
    -pe shmem 1 \
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

mkdir -p ${out_dir}
submit_knockout_job 0 0.01 "" "ptv,damaging_missense"
#submit_knockout_job 0 0.01 "" "ptv"
#submit_knockout_job 0 0.01 "" "synonymous"
#submit_knockout_job 0 0.01 "" "ptv,ptv_LC"

#submit_knockout_job "damaging_missense"
#submit_knockout_job "ptv"
#submit_knockout_job "synonymous"
#submit_knockout_job "ptv,ptv_LC"
#submit_knockout_job "ptv,ptv_LC,damaging_missense"


