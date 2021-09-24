#!/usr/bin/env bash
#
#
#$ -N _spa_test
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/step2_saige.log
#$ -e logs/step2_saige.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

# arguments
readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly in_vcf=${2?Error: Missing arg2 (in_vcf)}
readonly in_csi=${3?Error: Missing arg2 (in_csi)}
readonly in_gmat=${4?Error: Missing arg3 (in_gmat)} 
readonly in_var=${5?Error: Missing arg4 (in_var)} 
readonly out_prefix=${6?Error: Missing arg5 (path prefix for saige output)}
readonly chr=${SGE_TASK_ID}

# SAIGE paths
readonly threads=$(( ${NSLOTS}-1 ))
readonly createSparseGRM="/well/lindgren/flassen/software/dev/SAIGE/extdata/createSparseGRM.R"
readonly step1_fitNULLGLMM="/well/lindgren/flassen/software/dev/SAIGE/extdata/step1_fitNULLGLMM.R"
readonly step2_SPAtests="/well/lindgren/flassen/software/dev/SAIGE/extdata/step2_SPAtests.R"

# out prefix
readonly out="${out_prefix}_${chr}"

spa_test() {
   SECONDS=0
   print_update "starting SPA test for ${out}"
   set -x
   Rscript "${step2_SPAtests}"	\
	 --vcfFile=${in_vcf} \
	 --vcfFileIndex=${in_csi} \
	 --vcfField="DS" \
     --chrom=${chr} \
     --minMAF=0.00001 \
     --minMAC=1 \
     --GMMATmodelFile=${in_gmat} \
     --varianceRatioFile=${in_var} \
     --SAIGEOutputFile=${out} \
     --numLinesOutput=2 \
     --IsOutputAFinCaseCtrl=TRUE	\
	 --IsOutputNinCaseCtrl=TRUE	\
  	 --IsOutputHetHomCountsinCaseCtrl=TRUE
   set +x
   duration=${SECONDS}
   	print_update "Finished SPA test for ${out} at ${duration}"
}

if [ ! -f ${out}* ]; then
   set_up_RSAIGE
   spa_test
else
   print_update "${out} already exists. Skipping." | tee /dev/stderr
fi 


