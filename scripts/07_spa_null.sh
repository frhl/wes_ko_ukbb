#!/usr/bin/env bash
#
#$ -N submit_spa_null
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/submit_spa_null.log
#$ -e logs/submit_spa_null.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 3
#$ -V

module purge
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly plink_dir="data/saige/grm/input"
readonly grm_dir="data/saige/grm/input"
readonly covar_dir="data/phenotypes"
readonly pheno_dir="data/phenotypes"

readonly grm_mtx="${grm_dir}/211102_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"
readonly plink_file="${plink_dir}/211102_long_ukb_wes_200k_sparse_autosomes"
readonly covar_file="${covar_dir}/covars1.csv"
readonly covariates=$( cat ${covar_file} )
readonly pheno_file="${pheno_dir}/curated_phenotypes.tsv"
 
readonly spa_null_script="scripts/_spa_null.sh"

readonly queue="short.qe"
readonly nslots=3

fit_binary_traits() {
  local trait_type="binary"
  local out_dir="data/saige/output/combined/binary/step1"
  local pheno_list="${pheno_dir}/curated_phenotypes_binary_header.tsv"
  local phenotype=$( cut -f${SGE_TASK_ID} ${pheno_list} )
  local out_prefix="${out_dir}/ukb_wes_200k_${phenotype}"
  submit_spa_null

}

fit_cts_traits() {
  local trait_type="quantitative"
  local out_dir="data/saige/output/combined/cts/step1"
  local pheno_list="${pheno_dir}/curated_phenotypes_cts_header.tsv"
  local phenotype=$( cut -f${SGE_TASK_ID} ${pheno_list} )
  local out_prefix="${out_dir}/ukb_wes_200k_${phenotype}"
  submit_spa_null

}

submit_spa_null() {
  mkdir -p ${out_dir}
  if [ ! -f ${out_prefix}* ]; then
    set -x
    qsub -N "spa_${phenotype}" \
      -t "${SGE_TASK_ID}" \
      -q "${queue}" \
      -pe shmem ${nslots} \
      "${spa_null_script}" \
      "${plink_file}" \
      "${pheno_file}" \
      "${phenotype}" \
      "${covariates}" \
      "${trait_type}" \
      "${grm_mtx}" \
      "${grm_sam}" \
      "${out_prefix}"
    set +x
  else
    >&2 echo "${out_prefix} already exists. Skipping.."
  fi
}

# Fit null model for binary/cts traits
fit_binary_traits
fit_cts_traits





