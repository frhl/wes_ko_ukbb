#!/usr/bin/env bash
#
#
#$ -N _spa_test
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/spa_test.log
#$ -e logs/spa_test.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 22
#$ -q lindgren.qe
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

# arguments
readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly in_vcf=${2?Error: Missing arg2 (in_vcf)}
readonly in_csi=${2?Error: Missing arg2 (in_csi)}
readonly in_gmat=${2?Error: Missing arg3 (in_gmat)} 
readonly in_var=${2?Error: Missing arg4 (in_var)} 
readonly out_prefix=${8?Error: Missing arg5 (path prefix for saige output)}
readonly chr=${SGE_TASK_ID}

# SAIGE paths
readonly threads=$(( ${NSLOTS}-1 ))
readonly createSparseGRM="/well/lindgren/flassen/software/dev/SAIGE/extdata/createSparseGRM.R"
readonly step1_fitNULLGLMM="/well/lindgren/flassen/software/dev/SAIGE/extdata/step1_fitNULLGLMM.R"
readonly step2_SPAtests="/well/lindgren/flassen/software/dev/SAIGE/extdata/step2_SPAtests.R"

# out prefix
readonly out_prefix_chr="${out_prefix}_${chr}"

# run saige
set_up_RSAIGE
SECONDS=0
Rscript "${step1_SPAtests}"	\
	--vcfFile=${in_vcf} \
	--vcfFileIndex=${in_csi} \
	--vcfField="DS" \
    --chrom=${chr} \
    --minMAF=0.00001 \
    --minMAC=1 \
    --GMMATmodelFile=${in_gmat} \
    --varianceRatioFile=${in_var} \
    --SAIGEOutputFile=${out_prefix_chr} \
    --numLinesOutput=2 \
    --IsOutputAFinCaseCtrl=TRUE	\
	--IsOutputNinCaseCtrl=TRUE	\
	--IsOutputHetHomCountsinCaseCtrl=TRUE

duration=${SECONDS}
print_update "[${phenotype}-chr${chr}]: finished saige step 2 after ${duration}"


