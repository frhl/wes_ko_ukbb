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

readonly spark_dir="data/tmp/spark"
readonly bash_script="scripts/_knockouts.sh"

readonly in_dir="data/mt/annotated"
readonly out_dir="data/knockouts/alt"
readonly in_prefix="${in_dir}/ukb_eur_wes_200k_annot_chrCHR.mt"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k"
readonly out_type="vcf"

# Note: ~20 slots are needed for running chr1 
# when using aggr_method="collect" on short.qe
readonly tasks="21"
readonly queue="short.qa"
readonly nslots=1

readonly only_vcf=""
readonly aggr_method="fast" # either fasts or collect

# variant and sample parameters
readonly exclude="data/genes/220310_common_plofs_to_exclude.txt"
readonly maf_lb="0"
readonly maf_ub="5e-2"
readonly sex="both"

mkdir -p ${out_dir}

submit_knockout_job() 
{
  local annotation=${1}
  local input_type="mt"
  local input_path="${in_prefix}_chrCHR_maf${maf_lb}to${maf_ub}_${annotation}.mt"
  local prefix="${out_prefix}_chrCHR_maf${maf_lb}to${maf_ub}_${annotation}"
  set -x
  qsub -N "_ko_${csqs}" \
    -o "logs/_knockouts.log" \
    -e "logs/_knockouts.errors.log" \
    -t ${tasks} \
    -q "${queue}" \
    -pe shmem ${nslots} \
    "${bash_script}" \
    "${input_path}" \
    "${input_type}" \
    "${af_min}" \
    "${af_max}" \
    "${maf_lb}" \
    "${maf_ub}" \
    "${exclude}" \
    "${sex}" \
    "${annotation}" \
    "${only_vcf}" \
    "${aggr_method}" \
    "${prefix}" \
    "${out_type}" 
  set +x
}

#submit_knockout_job "pLoF_damaging_missense"
submit_knockout_job "pLoF"
#submit_knockout_job "0" "5e-2" "" "damaging_missense"
#submit_knockout_job "0" "5e-2" "" "synonymous"
#submit_knockout_job "0" "5e-2" "" "pLoF,LC,damaging_missense"

#submit_knockout_job 0 0.05 "" "pLoF"
#submit_knockout_job 0 0.05 "" "synonymous"
#§submit_knockout_job 0 0.05 "" "ptv,LC"
