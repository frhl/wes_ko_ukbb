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
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/submit_knockout.log
#$ -e logs/submit_knockout.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qf
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

# directories
readonly in_dir="data/mt"
readonly spark_dir="data/tmp/spark"
readonly out_dir="derived/knockouts/211206"

# input path
readonly knockout_script="scripts/_knockouts.sh"
readonly in_phased="${in_dir}/ukb_wes_200k_annotated_chrCHR.mt"
readonly in_unphased="${in_dir}/ukb_wes_200k_annotated_chrCHR_singletons.mt"
readonly in_phased_type="mt"
readonly in_unphased_type="mt"

# parameters required
readonly af_min=""
readonly af_max=""

# output path
readonly prefix="${out_dir}/ukb_wes_200k"

submit_knockout_job() 
{
  out_prefix="${prefix}_maf${1}_${2}_${3:+_${3}}chrCHR"
  qsub_name=$( echo ${4} | tr "," "_")
  set -x
  qsub -N "_ko_${qsub_name}" \
    -t 21 \
    -q "short.qa" \
    -pe shmem 2 \
    "${knockout_script}" \
    "${in_phased}" \
    "${in_phased_type}" \
    "${in_unphased}" \
    "${in_unphased_type}" \
    "${af_min}" \
    "${af_max}" \
    "${1}" \
    "${2}" \
    "${3}" \
    "${4}"\
    "${out_prefix}"
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


