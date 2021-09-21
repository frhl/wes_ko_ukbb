#!/usr/bin/env bash
#
#
#$ -N step2_saige
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/step2_saige.log
#$ -e logs/step2_saige.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qe
#$ -t 22

#set -o errexit
#set -o nounset

module purge
source utils/bash_utils.sh

SGE_TASK_ID=22

# directories
readonly out_dir="derived/saige/binary/step2"

# input path
readonly chr=${SGE_TASK_ID}
readonly covar_file="${covar_dir}/COVARS1.csv"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_test_crohns${chr}"
readonly out_prefix_ratio="${out_dir}/ukb_wes_200k_test_crohns_cate"

# SAIGE paths
readonly threads=$(( ${NSLOTS}-1 ))
readonly createSparseGRM="/well/lindgren/flassen/software/dev/SAIGE/extdata/createSparseGRM.R"
readonly step1_fitNULLGLMM="/well/lindgren/flassen/software/dev/SAIGE/extdata/step1_fitNULLGLMM.R"
readonly step2_SPAtests="/well/lindgren/flassen/software/dev/SAIGE/extdata/step2_SPAtests.R"

Rscript "${step1_SPAtests}"	\
	--vcfFile=./input/dosage_10markers.vcf.gz \
	--vcfFileIndex=./input/dosage_10markers.vcf.gz.tbi \
	--vcfField="DS" \
    --chrom=${chr} \
    --minMAF=0.0001 \
    --minMAC=1 \
    --GMMATmodelFile=./output/example_binary_includenonAutoforvarRatio.rda \
    --varianceRatioFile=./output/example_binary.varianceRatio.txt \
    --SAIGEOutputFile=./output/example_binary.SAIGE.vcf.genotype.txt_new \
    --numLinesOutput=2 \
    --IsOutputAFinCaseCtrl=TRUE	\
	--IsOutputNinCaseCtrl=TRUE	\
	--IsOutputHetHomCountsinCaseCtrl=TRUE	


