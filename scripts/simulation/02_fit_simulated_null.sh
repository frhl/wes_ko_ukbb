#!/usr/bin/env bash
#
#$ -N fit_simulated_null
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/fit_simulated_null.log
#$ -e logs/fit_simulated_null.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hge
#$ -t 1-10
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly plink_dir="data/saige/grm/input"
readonly grm_dir="data/saige/grm/input"
readonly covar_dir="data/phenotypes"
readonly pheno_dir="data/simulation/phenotypes"
readonly out_dir="data/simulation/saige/step1"

readonly grm_mtx="${grm_dir}/211102_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"

readonly plink_file="${plink_dir}/chunks/ukb_wes_200k_sparse_autosomes_mrg"
readonly covar_file="${covar_dir}/covars2.csv"
readonly covariates=$( cat ${covar_file} )

readonly step1_fitNULLGLMM="utils/saige/step1_fitNULLGLMM.R"


fit_cts() {
  local pheno_file=${1}
  local prefix=${2}
  local trait_type="quantitative"
  local inv_normalize="TRUE"
  local phenotype="y_cts_${SGE_TASK_ID}"
 fit_null
}

fit_bin() {
  local pheno_file=${1}
  local prefix=${2}
  local trait_type="binary"
  local inv_normalize="FALSE"
  local phenotype="y_bin_${SGE_TASK_ID}"
  fit_null 
}

fit_null() {
   mkdir -p ${out_dir}
   local out_prefix="${out_dir}/${prefix}_${phenotype}"
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
fit_cts "${pheno_dir}/ukb_eur_h2_0_0_pi_NA_NA_K_1e-1_chr21_phenos.tsv.gz" "ukb_eur_h2_0_0_pi_NA_NA_K_1e-1_chr21"
fit_cts "${pheno_dir}/ukb_eur_h2_0_3e-1_pi_NA_NA_K_1e-1_chr21_phenos.tsv.gz" "ukb_eur_h2_0_3e-1_pi_NA_NA_K_1e-1_chr21"
fit_cts "${pheno_dir}/ukb_eur_h2_2e-1_3e-1_pi_NA_NA_K_1e-1_chr21_phenos.tsv.gz" "ukb_eur_h2_2e-1_3e-1_pi_NA_NA_K_1e-1_chr21"
#fit_bin "${pheno_dir}/ukb_eur_h2_0_0_pi_NA_NA_K_1e-1_chr21_phenos.tsv.gz" "ukb_eur_h2_0_0_pi_NA_NA_K_1e-1_chr21"







