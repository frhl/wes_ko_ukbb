#!/usr/bin/env bash
#
# Perform a conditional analysis on pre-selected common markers.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=spa_cond_common
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/spa_cond_common.log
#SBATCH --error=logs/spa_cond_common.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-320
# --begin=now+2hour
# --dependency=afterok:13640293

set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly cluster=$( get_current_cluster)
readonly index=$( get_array_task_id )

readonly curwd="$(pwd)"
readonly vcf_dir="data/conditional/common/combined"
readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"
readonly rscript="scripts/_check_prs_ok.R"

readonly spa_script="scripts/conditional/common/_spa_cond_common.sh"
readonly merge_script="scripts/_spa_merge.sh"
readonly in_prefix="ukb_eur_wes_200k"

readonly cond_dir="data/conditional/common/combined"
readonly cond_file="${cond_dir}/ukb_eur_wes_200k_chrCHR_pLoF_damaging_missense_w_common_markers.txt"
readonly cond_cat="common" 

readonly grm_dir="data/saige/grm/input/dnanexus"
readonly grm_mtx="${grm_dir}/ukb_eur_200k_grm_fitted_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"
readonly plink_file="${grm_dir}/ukb_eur_200k_grm_grch38_rv_merged"

readonly markers_dir="data/conditional/common/markers"
readonly markers="${markers_dir}/common_conditional.markers"

submit_spa_binary_with_csqs()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/dec22_phenotypes_binary_200k_header.tsv"
  local phenotype=$( sed "${index}q;d" ${pheno_list} )
  submit_spa_with_csqs "${annotation}" "${phenotype}" "binary"
}

submit_spa_cts_with_csqs()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( sed "${index}q;d" ${pheno_list} )
  submit_spa_with_csqs "${annotation}" "${phenotype}" "cts"
}

submit_spa_with_csqs()
{
  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local trait=${3?Error: Missing arg3 (trait)}
  if [ ! -z ${phenotype} ]; then
    # only run phenotypes with actual common markers
    if [ "$(cat ${markers} | cut -f6 | uniq | grep -w ${phenotype} | wc -l )" -gt "0" ]; then

      local step1_dir="data/saige/output/${trait}/step1"
      local step2_dir="data/saige/output/${trait}/step2_common_new/min_mac${min_mac}"
      local in_vcf="${vcf_dir}/${in_prefix}_chrCHR_${annotation}_w_common.vcf.bgz"
      mkdir -p ${step2_dir}

      local in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
      local in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
      local out_prefix="${step2_dir}/${in_prefix}_chrCHR_${phenotype}_${annotation}"
      local out_mrg="${step2_dir}/${in_prefix}_${phenotype}_${annotation}.txt.gz"
      if [ "${use_prs}" -eq "1" ]; then
        set_up_rpy
        local in_gmat_prs="${step1_dir}/ukb_wes_200k_${phenotype}_chrCHR.rda"
        local in_var_prs="${step1_dir}/ukb_wes_200k_${phenotype}_chrCHR.varianceRatio.txt"
        local prs_ok=$(Rscript ${rscript} --phenotype ${phenotype})
        if [ -f "${in_gmat_prs/CHR/21}" ] & [ -f "${in_var_prs/CHR/21}" ] & [ "${prs_ok}" -eq "1" ]; then
           local in_gmat=${in_gmat_prs}
           local in_var=${in_var_prs}
           local out_prefix="${step2_dir}/${in_prefix}_chrCHR_${phenotype}_${annotation}_locoprs"
           local out_mrg="${step2_dir}/${in_prefix}_${phenotype}_${annotation}_locoprs.txt.gz"
         else
           >&2 echo "Using without PRS."
         fi
       fi

       if [ ! -f "${out_mrg}" ]; then
         local slurm_spa_name="spa_${phenotype}_${annotation}"
         local slurm_merge_name="_mrg_${phenotype}_${annotation}"
         submit_spa_job
         submit_merge_job
       else
         >&2 echo "Phenotype ${phenotype} with annotation ${annotation} already exists! Skipping.." 
       fi
    else
      >&2 echo "Phenotype ${phenotype} does not have conditional markers. Skipping.."
    fi
  else
    >&2 echo "No phenotype at index ${index}. Exiting.." 
  fi 
}


submit_spa_job() {
  mkdir -p ${step2_dir}
  local slurm_tasks="${tasks}"
  local slurm_jname="${slurm_spa_name}"
  local slurm_lname="logs/_spa_cond_common"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local slurm_nslots="${nslots}"
  readonly spa_jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname}.log" \
    --error="${slurm_lname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --array=${slurm_tasks} \
    --parsable \
    "${spa_script}" \
    "${phenotype}" \
    "${in_vcf}" \
    "${in_vcf}.csi" \
    "${in_gmat}" \
    "${in_var}" \
    "${grm_mtx}" \
    "${grm_sam}" \
    "${min_mac}" \
    "${out_prefix}" \
    "${cond_file}" \
    "${cond_cat}" )
}


submit_merge_job()
{
  local remove_by_chr="Y"
  local slurm_jname="${slurm_merge_name}"
  local slurm_lname="${out_prefix}"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local slurm_nslots="1"
  readonly merge_jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname}.log" \
    --error="${slurm_lname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --dependency="afterok:${spa_jid}" \
    --open-mode="append" \
    --parsable \
    "${merge_script}" \
    "${out_prefix}" \
    "${out_mrg}" \
    "${remove_by_chr}" )
}



# parameters
readonly use_prs="1"
readonly min_mac=4
readonly tasks=1-22
readonly queue="short"
readonly project="lindgren.prj"
readonly nslots=2



# cts traits
#submit_spa_cts_with_csqs "pLoF_damaging_missense"
submit_spa_binary_with_csqs "pLoF_damaging_missense"

#sleep 10
#submit_spa_cts_with_csqs "pLoF"
#submit_spa_binary_with_csqs "pLoF"

#sleep 10
#submit_spa_cts_with_csqs "damaging_missense"
#submit_spa_binary_with_csqs "damaging_missense"

#sleep 10
#submit_spa_cts_with_csqs "synonymous"
#submit_spa_binary_with_csqs "synonymous"





