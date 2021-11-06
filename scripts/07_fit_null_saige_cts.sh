#!/usr/bin/env bash
#
#$ -N saige_null_cts
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/saige_null_cts.log
#$ -e logs/saige_null_cts.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qe
#$ -t 1-40

module purge
source utils/bash_utils.sh
source utils/hail_utils.sh

# SAIGE paths
readonly threads=$(( ${NSLOTS}-1 ))
readonly createSparseGRM="/well/lindgren/flassen/software/dev/SAIGE/extdata/createSparseGRM.R"
readonly step1_fitNULLGLMM="/well/lindgren/flassen/software/dev/SAIGE/extdata/step1_fitNULLGLMM.R"
readonly step2_SPAtests="/well/lindgren/flassen/software/dev/SAIGE/extdata/step2_SPAtests.R"

# directories
readonly plink_dir="data/saige/grm/input"
readonly grm_dir="data/saige/grm/input"
readonly covar_dir="data/phenotypes"
readonly pheno_dir="data/phenotypes"
readonly out_dir="data/saige/output/combined/cts/step1"

# input path
readonly chr=${SGE_TASK_ID}
readonly grm_mtx="${grm_dir}/211102_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"
readonly plink_file="${plink_dir}/211102_long_ukb_wes_200k_sparse_autosomes"
readonly pheno_file="${pheno_dir}/UKBB_WES200k_filtered_cts_phenotypes.tsv.gz"
readonly covar_file="${covar_dir}/COVARS1_BMI.csv"
readonly pheno_list="${pheno_dir}/UKBB_WES200k_cts_phenotypes_header.txt"

# select covars and phenotype (0-42)
readonly index=${SGE_TASK_ID}
readonly covariates=$( cat ${covar_file} )
readonly phenotype=$( cut -f${index} ${pheno_list} )

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_adjBMI_${phenotype}"

# null model script
fit_null() {
   SECONDS=0
   print_update "Fitting null model for ${phenotype}, out: ${out_prefix}"
   Rscript "${step1_fitNULLGLMM}" \
     --plinkFile="${plink_file}" \
     --phenoFile="${pheno_file}" \
     --phenoCol="${phenotype}" \
     --covarColList=${covariates} \
     --sampleIDColinphenoFile="ID" \
     --traitType="quantitative" \
     --invNormalize=TRUE \
     --outputPrefix="${out_prefix}" \
     --outputPrefix_varRatio="${out_prefix}" \
     --IsOverwriteVarianceRatioFile=TRUE \
     --sparseGRMFile=${grm_mtx} \
     --sparseGRMSampleIDFile=${grm_sam}  \
     --nThreads=${threads} \
     --LOCO=FALSE\
     --skipModelFitting=FALSE \
     --IsSparseKin=TRUE      \
     --isCateVarianceRatio=FALSE # Only needed for SAIGE-GENE+ (geneset)	
   duration=${SECONDS}
   print_update "Finished fitting null model for ${phenotype}" "$duration"
}

if [ ! -f ${out_prefix}* ]; then
   set_up_RSAIGE
   fit_null
else
   print_update "Warning: Null model at ${out_prefix} already exists. Skipping."
fi 



