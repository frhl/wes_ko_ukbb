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

