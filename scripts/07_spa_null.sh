#!/usr/bin/env bash
#
#$ -N spa_null
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/spa_null.log
#$ -e logs/spa_null.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1-20
#$ -tc 1
#$ -V

# all binary: 1 - 71

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spa_null_script="scripts/_spa_null.sh"
readonly rscript="scripts/_spa_null.R"
readonly rscript_ldsc="scripts/_check_prs_p.R"

readonly plink_dir="data/saige/grm/input"
readonly grm_dir="data/saige/grm/input"
readonly covar_dir="data/phenotypes"
readonly pheno_dir="data/phenotypes"
readonly prs_dir="data/prs/scores"
readonly ldsc_dir="data/prs/ldsc"

readonly grm_mtx="${grm_dir}/211102_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"
#readonly plink_file="${plink_dir}/211102_long_ukb_wes_200k_sparse_autosomes"
readonly plink_file="${plink_dir}/chunks/ukb_wes_200k_sparse_autosomes_mrg"
readonly covar_file="${covar_dir}/covars1.csv"
readonly covariates=$( cat ${covar_file} )

readonly out_prefix="ukb_wes_200k"
readonly index=${SGE_TASK_ID}



fit_binary_traits() {
  local trait_type="binary"
  local inv_normalize="FALSE"
  local out_dir="data/saige/output/binary/step1"
  local pheno_list="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
  local phenotype=$( sed "${index}q;d" ${pheno_list} )
  local out="${out_dir}/${out_prefix}_${phenotype}"
  pheno_file="${pheno_dir}/filtered_covar_phenotypes_binary.tsv.gz"
  
  local out_pheno_prs="${out_dir}/${phenotype}_prs.txt.gz"
  local prs="${prs_dir}/${phenotype}_pgs_chrom.txt.gz"
  local ldsc="${ldsc_dir}/ldsc_${phenotype}.rds"
  submit_spa_null
}

fit_cts_traits() {
  local trait_type="quantitative"
  local inv_normalize="FALSE"
  local out_dir="data/saige/output/cts/step1"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype_in=$( sed "${index}q;d" ${pheno_list} )
  local phenotype="${phenotype_in}_int" # hack to add inverse norm transformed
  local out="${out_dir}/${out_prefix}_${phenotype}"
  pheno_file="${pheno_dir}/filtered_covar_phenotypes_cts.tsv.gz"

  local out_pheno_prs="${out_dir}/${phenotype}_int_prs.txt.gz"
  local prs="${prs_dir}/${phenotype}_int_pgs_chrom.txt.gz"
  local ldsc="${ldsc_dir}/ldsc_${phenotype}_int.rds"
  submit_spa_null
}

set_up_prs() {
  # assuming that this function is called downstream of 
  # fit_cts_traits/fit_binary_traits
  set_up_rpy
  prs_ok=0
  if [[ -f "${prs}"  && "${use_prs}" -eq "1" ]]; then
    echo "Note: Checking LDSC h2 estimates at ${ldsc}."
    if [ -f "${ldsc}" ]; then
      local h2_pass_qc=$(Rscript ${rscript_ldsc} --ldsc ${ldsc})
      echo "Note: PRS for ${phenotype} pass QC: ${h2_pass_qc}"
      if [ "${h2_pass_qc}" -eq "1" ]; then
        if [ ! -f "${out_pheno_prs}" ]; then
          Rscript ${rscript} \
            --phenotype ${phenotype} \
            --covariates ${covariates} \
            --phenofile ${pheno_file} \
            --prsfile ${prs} \
            --outfile ${out_pheno_prs}
        fi
        # if the pheno file has already been created
        # ensure that chromosomal tasks are set up
        prs_ok=1
        pheno_file=${out_pheno_prs}
        tasks=1-22
      fi
    fi
  fi
}


submit_spa_null() {
  mkdir -p ${out_dir}
  tasks=${SGE_TASK_ID}
  set_up_prs
  if [[ "${prs_ok}" -eq "0"  && "${use_prs}" -eq "1" ]]; then
    >&2 echo "Note: PRS could not be started for ${phenotype}. Skipping."
  else
    if [ ! -z ${phenotype} ]; then
      if [ ! -f "${out_prefix}.rda" ]; then
        set -x
        qsub -N "_null_${phenotype}" \
         -t "${tasks}" \
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
          "${inv_normalize}" \
          "${use_prs}" \
          "${out}"
        set +x
      else
        >&2 echo "${out_prefix} already exists. Skipping.."
      fi
    else
      >&2 echo "No phenotype at index ${SGE_TASK_ID}. Exiting.." 
    fi
  fi
}

# Parameters
readonly use_prs=1
readonly nslots=2
readonly queue="short.qe"

# Fit null model for binary/cts traits
fit_cts_traits
fit_binary_traits





