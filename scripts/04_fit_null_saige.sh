#!/usr/bin/env bash
#
#$ -N step1_saige
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/step1_saige.log
#$ -e logs/step1_saige.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 8
#$ -q short.qe
#$ -t 1-42


module purge
source utils/bash_utils.sh
source utils/hail_utils.sh

# SAIGE paths
readonly threads=$(( ${NSLOTS}-1 ))
readonly createSparseGRM="/well/lindgren/flassen/software/dev/SAIGE/extdata/createSparseGRM.R"
readonly step1_fitNULLGLMM="/well/lindgren/flassen/software/dev/SAIGE/extdata/step1_fitNULLGLMM.R"
readonly step2_SPAtests="/well/lindgren/flassen/software/dev/SAIGE/extdata/step2_SPAtests.R"

# directories
readonly grm_dir="data/saige/grm/input"
readonly covar_dir="data/saige/input"
readonly pheno_dir="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/200k"
readonly plink_dir="data/saige/grm/input"
readonly out_dir="derived/saige/binary/step1"

# input path
readonly chr=${SGE_TASK_ID}
readonly grm_mtx="${grm_dir}/ukb_eur_wes_sparse_markers_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_dir}/ukb_eur_wes_sparse_markers_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
readonly plink_file="${plink_dir}/ukb_eur_wes_sparse_markers_autosomes"
readonly pheno_file="${pheno_dir}/UKBB_WES200k_filtered_binary_phenotypes.tsv"
readonly covar_file="${covar_dir}/COVARS1.csv"
readonly pyscript="utils/subscript/extract_phenos_from_header.py"

# select covars and phenotype (0-42)
set_up_pythonpath
covariates=$( cat ${covar_file} )
phenotype=$( python ${pyscript} \
    --input ${pheno_file} \
    --index ${SGE_TASK_ID})

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_${phenotype}"

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
     --traitType="binary" \
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



