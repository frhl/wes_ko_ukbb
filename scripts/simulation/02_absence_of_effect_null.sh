#!/usr/bin/env bash
#
#$ -N absence_of_effect_null
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/absence_of_effect_null.log
#$ -e logs/absence_of_effect_null.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hge
#$ -t 2
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly plink_dir="data/saige/grm/input"
readonly grm_dir="data/saige/grm/input"
readonly covar_dir="data/phenotypes"
readonly pheno_dir="data/simulation/absence_of_effect"
readonly out_dir="data/simulation/saige/step1"

readonly grm_mtx="${grm_dir}/211102_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"

readonly plink_file="${plink_dir}/chunks/ukb_wes_200k_sparse_autosomes_mrg"
readonly covar_file="${covar_dir}/covars1.csv"
readonly covariates=$( cat ${covar_file} )

readonly chr="21"
readonly pheno_file="${pheno_dir}/ukb_eur_h2_0_pi_None_chr${chr}_phenotype.tsv.gz"
readonly step1_fitNULLGLMM="utils/saige/step1_fitNULLGLMM.R"


fit_cts() {
  local trait_type="quantitative"
  local inv_normalize="TRUE"
  local phenotype="cts${SGE_TASK_ID}"
 fit_null 
}

fit_bin() {
  local trait_type="binary"
  local inv_normalize="FALSE"
  local phenotype="bin${SGE_TASK_ID}"
  fit_null 
}

fit_null() {
   mkdir -p ${out_dir}
   local out_prefix="${out_dir}/ukb_eur_h2_0_pi_None_${phenotype}_chr${chr}"
   if [ ! -f "${out_prefix}.rda" ]; then
     SECONDS=0
     Rscript "${step1_fitNULLGLMM}" \
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
 }


# Running null model
set_up_RSAIGE
fit_cts



