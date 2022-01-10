#!/usr/bin/env bash
#
#
#$ -N _spa_test
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/spa_test.log
#$ -e logs/spa_test.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly in_vcf=${2?Error: Missing arg2 (in_vcf)}
readonly in_csi=${3?Error: Missing arg3 (in_csi)}
readonly in_gmat=${4?Error: Missing arg4 (in_gmat)} 
readonly in_var=${5?Error: Missing arg5 (in_var)} 
readonly out_prefix=${6?Error: Missing arg6 (path prefix for saige output)}
readonly in_markers=${7} # optional conditioning markers
readonly chr=${SGE_TASK_ID}

readonly vcf=$(echo ${in_vcf} | sed -e "s/CHR/${chr}/g")
readonly csi=$(echo ${in_csi} | sed -e "s/CHR/${chr}/g")
readonly markers=$(cat $(echo ${in_markers} | sed -e "s/CHR/${chr}/g"))

readonly threads=$(( ${NSLOTS}-1 ))
readonly createSparseGRM="/well/lindgren/flassen/software/dev/SAIGE/extdata/createSparseGRM.R"
readonly step1_fitNULLGLMM="/well/lindgren/flassen/software/dev/SAIGE/extdata/step1_fitNULLGLMM.R"
readonly step2_SPAtests="/well/lindgren/flassen/software/dev/SAIGE/extdata/step2_SPAtests.R"

readonly out="${out_prefix}_chr${chr}"

spa_test() {
   SECONDS=0
   print_update "SAIGE-SPA: Writing '${out}' using '${vcf}'"
   set -x
   Rscript "${step2_SPAtests}"	\
     --vcfFile=${vcf} \
     --vcfFileIndex=${csi} \
     --vcfField="DS" \
     --chrom="chr${chr}" \
     --minMAF=0.0000001 \
     --minMAC=1 \
     --GMMATmodelFile=${in_gmat} \
     --varianceRatioFile=${in_var} \
     --SAIGEOutputFile=${out} \
     --numLinesOutput=2 \
     --IsOutputAFinCaseCtrl=TRUE \
     --IsOutputNinCaseCtrl=TRUE \
     --IsOutputHetHomCountsinCaseCtrl=TRUE \
     --LOCO=FALSE\
     ${markers:+--condition "$markers"}
   set +x
   duration=${SECONDS}
   	print_update "Finished SPA test for ${out} at ${duration}"
}

if [ ! -f ${out}* ]; then
   set_up_RSAIGE
   rm ${out}*
   spa_test
else
   print_update "${out} already exists. Skipping." | tee /dev/stderr
fi 


