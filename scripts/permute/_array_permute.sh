#!/usr/bin/env bash
#
#
#$ -N _array_permute
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_array_permute.log
#$ -e logs/_array_permute.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly in_vcf=${1?Error: Missing arg1 (phenotype)}
readonly in_tsv=${2?Error: Missing arg2 (in_vcf)}
readonly in_spa=${3?Error: Missing arg3 (in_csi)}
readonly out_prefix=${7?Error: Missing arg6 (path prefix for saige output)}
readonly chr=${SGE_TASK_ID}

readonly vcf=$(echo ${in_vcf} | sed -e "s/CHR/${chr}/g")
readonly out=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")

readonly threads=$(( ${NSLOTS}-1 ))
readonly step2_SPAtests="utils/saige/step2_SPAtests.R"


calc_permute_chunk() {



}



submit_permute_chunk() {


}

if [ ! -f ${out} ]; then
   set_up_RSAIGE
   spa_test
else
  >&2 echo "${out} already exists. Skipping.."
fi 


