#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly cluster=$( get_current_cluster )
readonly index=$( get_array_task_id )
readonly threads=$( get_threads )

readonly plink_file=${1?Error: Missing arg1 (plink_file)}
readonly pheno_file=${2?Error: Missing arg2 (pheno_file)}
readonly phenotype=${3?Error: Missing arg3 (phenotype)}
readonly in_covariates=${4?Error: Missing arg4 (covariates)}
readonly trait_type=${5?Error: Missing arg5 (trait_type)}
readonly grm_mtx=${6?Error: Missing arg6 (grm_mtx)}
readonly grm_sam=${7?Error: Missing arg7 (grm_sam)}
readonly inv_normalize=${8?Error: Missing arg8(inv_normalize)}
readonly use_loco_prs=${9?Error: Missing arg8(inv_normalize)}
readonly out_prefix=${10?Error: Missing arg10 (out_prefix)}

readonly step1_fitNULLGLMM="utils/saige/step1_fitNULLGLMM.R"


fit_null() {
   if [ ! -f "${real_out_prefix}.rda" ]; then
     SECONDS=0
     set -x
     Rscript "${step1_fitNULLGLMM}" \
       --plinkFile="${plink_file}" \
       --phenoFile="${pheno_file}" \
       --phenoCol="${phenotype}" \
       --covarColList=${covariates} \
       --qCovarColList="sex,ukbb.centre"\
       --sampleIDColinphenoFile="eid" \
       --traitType="${trait_type}" \
       --invNormalize="${inv_normalize}" \
       --outputPrefix="${real_out_prefix}" \
       --outputPrefix_varRatio="${real_out_prefix}" \
       --IsOverwriteVarianceRatioFile=TRUE \
       --isCovariateOffset=FALSE \
       --sparseGRMFile=${grm_mtx} \
       --sparseGRMSampleIDFile=${grm_sam}  \
       --nThreads=${threads} \
       --LOCO=FALSE \
       --useSparseGRMtoFitNULL=TRUE \
       --isCateVarianceRatio=TRUE
   set +x
   else
     >&2 echo "Warning: Null model at ${real_out_prefix} already exists. Skipping!"
   fi
 }

# generate string of parmaeters used for off-chromosome PRS
get_loco_seq(){
  local target=${1}
  local sequence=$( for i in `seq 1 22`; do echo "chr$i"; done )
  if [ "$( echo ${sequence} | grep ${target} | wc -l)" -eq "1" ];then 
    local loco=$( echo ${sequence/"${target}"/""} | tr " " "," )
    echo ${loco}
  else
    raise_error "${target} is not a valid chromosome!"
  fi
}


# set up LOCO PRS conditioning
if [ "${use_loco_prs}" -eq "1" ]; then
  readonly chr="chr${index}"
  readonly loco_chr="loco_${chr}"
  readonly covariates="${in_covariates},${loco_chr}"
  readonly real_out_prefix="${out_prefix}_${chr}"
else
  readonly covariates="${in_covariates}"
  readonly real_out_prefix="${out_prefix}"
fi

# run model
set_up_RSAIGE
fit_null


