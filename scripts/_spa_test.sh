#!/usr/bin/env bash
#
#
#$ -N _spa_test
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_spa_test.log
#$ -e logs/_spa_test.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly in_vcf=${2?Error: Missing arg2 (in_vcf)}
readonly in_csi=${3?Error: Missing arg3 (in_csi)}
readonly in_gmat=${4?Error: Missing arg4 (in_gmat)} 
readonly in_var=${5?Error: Missing arg5 (in_var)} 
readonly min_mac=${6?Error: Missing arg6 (min_mac)} 
readonly out_prefix=${7?Error: Missing arg6 (path prefix for saige output)}
readonly in_markers="${8}" # optional conditioning markers
readonly chr=${SGE_TASK_ID}

readonly vcf=$(echo ${in_vcf} | sed -e "s/CHR/${chr}/g")
readonly csi=$(echo ${in_csi} | sed -e "s/CHR/${chr}/g")
readonly markers=$(cat $(echo ${in_markers} | sed -e "s/CHR/${chr}/g"))
readonly out=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")

readonly threads=$(( ${NSLOTS}-1 ))
readonly step2_SPAtests="utils/saige/step2_SPAtests.R"

spa_test() {
   SECONDS=0
   Rscript "${step2_SPAtests}"	\
     --vcfFile=${vcf} \
     --vcfFileIndex=${csi} \
     --vcfField="DS" \
     --chrom="chr${chr}" \
     --minMAF=0.0000001 \
     --minMAC=${min_mac} \
     --GMMATmodelFile=${in_gmat} \
     --varianceRatioFile=${in_var} \
     --SAIGEOutputFile=${out} \
     --LOCO=FALSE\
     ${markers:+--condition "$markers"} \
     && print_update "Finished saddle-point approximation for chr${chr}" ${SECONDS} \
     || raise_error "Saddle-point approximation for chr${chr} failed"
}

if [ ! -f ${out} ]; then
   set_up_RSAIGE
   spa_test
else
  >&2 echo "${out} already exists. Skipping.."
fi 


