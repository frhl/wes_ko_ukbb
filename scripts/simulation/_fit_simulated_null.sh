#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly rscript=${1?Error: Missing arg1 (phenotype)}
readonly plink_file=${2?Error: Missing arg2 (in_vcf)}
readonly pheno_file=${3?Error: Missing arg3 ()}
readonly in_phenotype=${4?Error: Missing arg3 ()}
readonly covariates=${5?Error: Missing arg3 ()}
readonly trait_type=${6?Error: Missing arg3 ()}
readonly inv_normalize=${7?Error: Missing arg3 ()}
readonly in_prefix=${8?Error: Missing arg3 ()}
readonly grm_mtx=${9?Error: Missing arg3 ()}
readonly grm_sam=${10?Error: Missing arg3 ()}

readonly array_idx=$( get_array_task_id )
readonly phenotype="${in_phenotype}_${array_idx}"
readonly out_prefix="${in_prefix}_${array_idx}"

if [ ! -f "${out_prefix}.rda" ]; then
 SECONDS=0
 set_up_RSAIGE
 set -x
 Rscript "${rscript}" \
   --plinkFile="${plink_file}" \
   --phenoFile="${pheno_file}" \
   --phenoCol="${phenotype}" \
   --covarColList=${covariates} \
   --sampleIDColinphenoFile="eid" \
   --traitType="${trait_type}" \
   --invNormalize="${inv_normalize}" \
   --outputPrefix="${out_prefix}" \
   --outputPrefix_varRatio="${out_prefix}" \
   --IsOverwriteVarianceRatioFile=TRUE \
   --isCovariateOffset=FALSE \
   --sparseGRMFile=${grm_mtx} \
   --sparseGRMSampleIDFile=${grm_sam}  \
   --nThreads=1 \
   --LOCO=FALSE \
   --useSparseGRMtoFitNULL=TRUE \
   --isCateVarianceRatio=FALSE \
   && print_update "Finished running SAIGE NULL model for ${phenotype}" ${SECONDS} \
   || raise_error "SAIGE NULL model failed for ${phenotype}"
else
 >&2 echo "Warning: Null model at ${out_prefix} already exists. Skipping!"
fi


