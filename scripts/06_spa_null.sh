#!/usr/bin/env bash
#
# @description generate saige null models (with and without PRS)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=spa_null
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/spa_null.log
#SBATCH --error=logs/spa_null.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-200


# all binary: 1 - 71

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly curwd=$(pwd)
readonly spa_null_script="scripts/_spa_null.sh"
readonly rscript="scripts/_spa_null.R"
readonly rscript_ldsc="scripts/_check_prs_p.R"

readonly plink_dir="data/saige/grm/input"
readonly grm_dir="data/saige/grm/input/dnanexus"
readonly covar_dir="data/phenotypes"
readonly pheno_dir="data/phenotypes"
readonly prs_dir="data/prs/scores"
readonly ldsc_dir="data/prs/ldsc"

readonly grm_mtx="${grm_dir}/ukb_eur_200k_grm_fitted_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"
readonly plink_file="${grm_dir}/ukb_eur_200k_grm_grch38_rv_merged"
readonly covar_file="${covar_dir}/covars1.csv"
readonly covariates=$( cat ${covar_file} )

readonly out_prefix="ukb_wes_200k"
readonly index=${SLURM_ARRAY_TASK_ID}



fit_binary_traits() {
  local trait_type="binary"
  local inv_normalize="FALSE"
  local out_dir="data/saige/output/binary/step1"
  local pheno_list="${pheno_dir}/spiros_brava_phenotypes_binary_200k_header.tsv"
  local phenotype=$( sed "${index}q;d" ${pheno_list} )
  local out="${out_dir}/${out_prefix}_${phenotype}"
  pheno_file="${pheno_dir}/spiros_brava_phenotypes_binary_200k.tsv"
  
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
  tasks=${index}
  set_up_prs
  if [[ "${prs_ok}" -eq "0"  && "${use_prs}" -eq "1" ]]; then
    >&2 echo "Note: PRS could not be started for ${phenotype}. Skipping."
  else
    if [ ! -z ${phenotype} ]; then
      if [ ! -f "${out_prefix}.rda" ]; then
        local slurm_tasks="${tasks}"
        local slurm_jname="_null_${phenotype}"
        local slurm_lname="logs/_spa_null"
        local slurm_project="${project}"
        local slurm_queue="${queue}"
        local slurm_nslots="${nslots}"
        readonly spa_null_jid=$( sbatch \
          --account="${slurm_project}" \
          --job-name="${slurm_jname}" \
          --output="${slurm_lname}.log" \
          --error="${slurm_lname}.errors.log" \
          --chdir="$(pwd)" \
          --partition="${slurm_queue}" \
          --cpus-per-task="${slurm_nslots}" \
          --array=${slurm_tasks} \
          --parsable \
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
          "${out}" )
      else
        >&2 echo "${out_prefix} already exists. Skipping.."
      fi
    else
      >&2 echo "No phenotype at index ${index}. Exiting.." 
    fi
  fi
}

# Parameters
readonly use_prs=0
readonly nslots=2
readonly queue="short"
readonly project="lindgren.prj"

# Fit null model for binary/cts traits
#fit_cts_traits
fit_binary_traits





