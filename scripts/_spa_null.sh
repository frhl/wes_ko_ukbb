#!/usr/bin/env bash
#
#$ -N _spa_null
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_spa_null.log
#$ -e logs/_spa_null.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly plink_file=${1?Error: Missing arg1 (phenotype)}
readonly pheno_file=${2?Error: Missing arg2 (in_vcf)}
readonly phenotype=${3?Error: Missing arg3 (in_csi)}
readonly covariates=${4?Error: Missing arg4 (in_gmat)} 
readonly trait_type=${5?Error: Missing arg5 (in_var)} 
readonly grm_mtx=${6?Error: Missing arg6 (path prefix for saige output)}
readonly grm_sam=${7?Error: Missing arg7 (Optional file with conditioning markers)}
readonly out_prefix=${8?Error: Missing arg7 (Optional file with conditioning markers)}

readonly createSparseGRM="/well/lindgren/flassen/software/dev/SAIGE/extdata/createSparseGRM.R"
readonly step1_fitNULLGLMM="/well/lindgren/flassen/software/dev/SAIGE/extdata/step1_fitNULLGLMM.R"
readonly threads=$(( ${NSLOTS}-1 ))

fit_null() {
   if [ ! -f "${out_prefix}.rda" ]; then
     SECONDS=0
     Rscript "${step1_fitNULLGLMM}" \
       --plinkFile="${plink_file}" \
       --phenoFile="${pheno_file}" \
       --phenoCol="${phenotype}" \
       --covarColList=${covariates} \
       --sampleIDColinphenoFile="eid" \
       --traitType="${trait_type}" \
       --invNormalize=TRUE \
       --outputPrefix="${out_prefix}" \
       --outputPrefix_varRatio="${out_prefix}" \
       --IsOverwriteVarianceRatioFile=TRUE \
       --sparseGRMFile=${grm_mtx} \
       --sparseGRMSampleIDFile=${grm_sam}  \
       --nThreads=${threads} \
       --LOCO=FALSE \
       --skipModelFitting=FALSE \
       --IsSparseKin=TRUE \
       --isCateVarianceRatio=FALSE \
       && print_update "Finished running SAIGE NULL model for ${phenotype}" ${SECONDS} \
       || raise_error "SAIGE NULL model failed for ${phenptype}"
       #--isCateVarianceRatio=FALSE  Only needed for SAIGE-GENE+ (geneset)
   else
     >&2 echo "Warning: Null model at ${out_prefix} already exists. Skipping!"
   fi
 }


# run NULL model
set_up_RSAIGE
fit_null


