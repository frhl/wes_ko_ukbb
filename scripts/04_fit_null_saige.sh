#!/usr/bin/env bash
#
#
#$ -N step1_saige
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/step1_saige.log
#$ -e logs/step1_saige.errors.log
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
readonly grm_dir="data/saige/grm/input"
readonly covar_dir="data/saige/input"
readonly pheno_dir="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/200k"
readonly plink_dir="data/saige/grm/input"
readonly out_dir="derived/saige/binary"

# input path
readonly chr=${SGE_TASK_ID}
readonly grm_mtx="${grm_dir}/ukb_imp_eur_chr1_22_sparse_markers_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_dir}/ukb_imp_eur_chr1_22_sparse_markers_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
readonly plink_file="${plink_dir}/ukb_imp_eur_chr1_22_sparse_markers"
readonly pheno_file="${pheno_dir}/UKBB_WES200k_binary_phenotypes.tsv"
readonly covar_file="${covar_dir}/COVARS1.csv"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_test_crohns${chr}"
readonly out_prefix_ratio="${out_dir}/ukb_wes_200k_test_crohns_cate"

# SAIGE paths
readonly threads=$(( ${NSLOTS}-1 ))
readonly createSparseGRM="/well/lindgren/flassen/software/dev/SAIGE/extdata/createSparseGRM.R"
readonly step1_fitNULLGLMM="/well/lindgren/flassen/software/dev/SAIGE/extdata/step1_fitNULLGLMM.R"
readonly step2_SPAtests="/well/lindgren/flassen/software/dev/SAIGE/extdata/step2_SPAtests.R"

covars=$( cat ${covar_file} )
set_up_RSAIGE
Rscript "${step1_fitNULLGLMM}" \
	--plinkFile="${plink_file}" \
    --phenoFile="${pheno_file}" \
    --phenoCol="Crohns_disease" \
    --covarColList=${covars} \
    --sampleIDColinphenoFile="ID" \
    --traitType="binary" \
    --invNormalize=TRUE \
    --outputPrefix="${out_prefix}" \
	--outputPrefix_varRatio="${out_prefix}" \
	--IsOverwriteVarianceRatioFile=TRUE \
    --sparseGRMFile=${grm_mtx} \
    --sparseGRMSampleIDFile=${grm_sam}  \
    --nThreads=1 \
    --LOCO=FALSE\
	--skipModelFitting=FALSE \
    --IsSparseKin=TRUE      \
    --isCateVarianceRatio=FALSE # Only needed for SAIGE-GENE+ (geneset)	

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


