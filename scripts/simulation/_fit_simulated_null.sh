#!/usr/bin/env bash
#
#$ -N _fit_simulated_null
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_fit_simulated_null.log
#$ -e logs/_fit_simulated_null.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hge
#$ -t 1
#$ -V

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

readonly phenotype="${in_phenotype}_${SGE_TASK_ID}"
readonly out_prefix="${in_prefix}_${SGE_TASK_ID}"

if [ ! -f "${out_prefix}.rda" ]; then
 SECONDS=0
 set_up_RSAIGE
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


