#!/usr/bin/env bash
#
#$ -N _spa_null
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/spa_null.log
#$ -e logs/spa_null.errors.log
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

# null model script
fit_null() {
   >&2 echo "${out_prefix}"
   if [ ! -f "${out_prefix}*" ]; then
     SECONDS=0
     set -x
     Rscript "${step1_fitNULLGLMM}" \
       --plinkFile="${plink_file}" \
       --phenoFile="${pheno_file}" \
       --phenoCol="${phenotype}" \
       --covarColList=${covariates} \
       --sampleIDColinphenoFile="ID" \
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
       --isCateVarianceRatio=FALSE # Only needed for SAIGE-GENE+ (geneset)
     set +x
     duration=${SECONDS}
     print_update "Finished fitting null model for ${phenotype}" "$duration"
   else
     print_update "Warning: Null model at ${out_prefix} already exists. Skipping."
   fi
 }


# run NULL model
set_up_RSAIGE
fit_null


