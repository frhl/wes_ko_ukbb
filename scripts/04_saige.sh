#!/usr/bin/env bash
#
#
#$ -N saige_step1
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/saige_step1.log
#$ -e logs/saige_step1.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qe
#$ -t 22

#set -o errexit
#set -o nounset

module purge
source utils/bash_utils.sh

SGE_TASK_ID=22

# directories
readonly grm_dir="data/saige/grm/input"
readonly covar_dir="data/saige"
readonly pheno_dir="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/200k"
readonly plink_dir="data/saige/grm/input"
readonly out_dir="derived/saige/binary"

# input path
readonly chr=${SGE_TASK_ID}
readonly grm_mtx="${grm_dir}/ukb_imp_eur_chr1_22_sparse_markers_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_dir}/ukb_imp_eur_chr1_22_sparse_markers_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
readonly plink_file="${plink_dir}/ukb_imp_eur_chr1_22_sparse_markers"
readonly pheno_file="${plink_dir}/UKBB_WES200k_binary_phenotypes.tsv"
readonly covar_file="${covar_dir}/COVAR_FILE_HERE"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_chr${chr}"

# SAIGE paths
readonly threads=$(( ${NSLOTS}-1 ))
readonly createSparseGRM="/well/lindgren/flassen/software/dev/SAIGE/extdata/createSparseGRM.R"
readonly step1_fitNULLGLMM="/well/lindgren/flassen/software/dev/SAIGE/step1_fitNULLGLMM.R"
readonly step2_SPAtests="/well/lindgren/flassen/software/dev/SAIGE/extdata/step2_SPAtests.R"


Rscript "${step1_fitNULLGLMM.R}" \
	--plinkFile="${plink_file}" \
    --phenoFile="${pheno_file}" \
    --phenoCol="Phenotype_here?" \
    --covarColList=${covar_file} \
    --sampleIDColinphenoFile="ID" \
    --traitType="binary" \
    --invNormalize=TRUE \
    --outputPrefix=./output/example_quantitative \
	--outputPrefix_varRatio=./output/example_quantitative_cate	\
	--sparseGRMFile=${grm_mtx} \
    --sparseGRMSampleIDFile=${grm_sam}  \
    --nThreads=4 \
    --LOCO=FALSE\
	--skipModelFitting=FALSE \
    --IsSparseKin=TRUE      \
    --isCateVarianceRatio=TRUE	

#Rscript step2_SPAtests.R	\
#	--vcfFile=./input/dosage_10markers.vcf.gz \
#	--vcfFileIndex=./input/dosage_10markers.vcf.gz.tbi \
#	--vcfField=DS \
#        --chrom=1 \
#        --minMAF=0.0001 \
#        --minMAC=1 \
#        --GMMATmodelFile=./output/example_binary_includenonAutoforvarRatio.rda \
#        --varianceRatioFile=./output/example_binary.varianceRatio.txt \
#        --SAIGEOutputFile=./output/example_binary.SAIGE.vcf.genotype.txt_new \
#        --numLinesOutput=2 \
#        --IsOutputAFinCaseCtrl=TRUE	\
#	--IsOutputNinCaseCtrl=TRUE	\
#	--IsOutputHetHomCountsinCaseCtrl=TRUE	


