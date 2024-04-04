#!/usr/bin/env bash
#
# condition on common markers near knockout genes that are significant 
# in primary analysis using an iterative approach. Note, that this script
# calls a sub script that should be in the long queue (probably due to the fact
# that it takes >24h to condition our signal in HLA region)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=spa_iter_common
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/spa_iter_common.log
#SBATCH --error=logs/spa_iter_common.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-350


set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly curwd=$(pwd)
readonly bash_script="scripts/conditional/common/_spa_iter_common.sh"
readonly rscript="scripts/_check_prs_ok.R"

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )

# parameters
readonly min_maf=0.01
readonly min_mac=4
readonly max_iter=25
readonly P_cutoff="5e-6"

# directories and paths
readonly pheno_dir="data/phenotypes"
readonly interval_dir="data/conditional/common/intervals/2024/min_mac${min_mac}"
readonly out_dir="data/conditional/common/spa_iter_new_2024"

readonly grm_dir="data/saige/grm/input/dnanexus"
readonly grm_mtx="${grm_dir}/ukb_eur_200k_grm_fitted_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"

readonly in_prefix="ukb_eur_wes_200k"

# use PRS and give error if not available
#readonly force_prs="0"
# use PRS of available
readonly use_prs="1"

submit_binary_analysis()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/dec22_phenotypes_binary_200k_header.tsv"
  local phenotype=$( sed "${task_id}q;d" ${pheno_list} )
  submit_cond_spa "${annotation}" "${phenotype}" "binary"
}

submit_cts_analysis()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( sed "${task_id}q;d" ${pheno_list} )
  submit_cond_spa "${annotation}" "${phenotype}" "cts"
}

submit_cond_spa()
{
  set_up_rpy
  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local trait=${3?Error: Missing arg3 (trait)}
  local out_prefix="${out_dir}/${in_prefix}_${phenotype}_${annotation}_cond"
  local step1_dir="data/saige/output/${trait}/step1"
  local step2_dir="data/saige/output/${trait}/step2/minmac${min_mac}"
  local in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
  local in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
  local prs_ok=$(Rscript ${rscript} --phenotype ${phenotype})
  if [ "${use_prs}" -eq "1" ] & [ "${prs_ok}" -eq "1" ]; then
     local in_gmat_prs="${step1_dir}/ukb_wes_200k_${phenotype}_chrCHR.rda"
     local in_var_prs="${step1_dir}/ukb_wes_200k_${phenotype}_chrCHR.varianceRatio.txt"
     if [ -f "${in_gmat_prs/CHR/21}" ] & [ -f "${in_var_prs/CHR/21}" ]; then
        local in_gmat=${in_gmat_prs}
        local in_var=${in_var_prs}
        local out_prefix="${out_prefix}_locoprs"
     else
        >&2 echo "Saige NULL (PRS) ${in_gmat_prs}/${in_var_prs} does not exist. Using without PRS."
     fi
  fi 
  

  
  mkdir -p ${out_dir}
  local intervals="${interval_dir}/${in_prefix}_${phenotype}_${annotation}_intervals.txt"
  local path_min_maf="${interval_dir}/${in_prefix}_${phenotype}_${annotation}_min_maf.tsv"
  local new_min_maf=$( sed "1q;d" ${path_min_maf} | cut -f4) 
  if [ -f "${intervals}" ] & [ -f "${path_min_maf}" ]; then 
    local n_vcfs=$( cat ${intervals} | wc -l )
    local tasks="1-${n_vcfs}"
    readonly slurm_jname="_cond_${phenotype}"
    readonly slurm_lname="${out_prefix}"
    readonly slurm_project="lindgren.prj"
    readonly slurm_queue="long"
    readonly sge_queue="long.qc"
    readonly slurm_shmem="1"
    if [ "${cluster}" = "slurm" ]; then
      sbatch \
        --account="${slurm_project}" \
        --job-name="${slurm_jname}" \
        --output="${slurm_lname}.log" \
        --error="${slurm_lname}.errors.log" \
        --chdir="${curwd}" \
        --partition="${slurm_queue}" \
        --cpus-per-task="${slurm_shmem}" \
        --array=${tasks} \
        --parsable \
        "${bash_script}" \
        "${in_gmat}" \
        "${in_var}" \
        "${intervals}" \
        "${out_prefix}" \
        "${P_cutoff}" \
        "${max_iter}" \
        "${min_mac}" \
        "${grm_mtx}" \
        "${grm_sam}" \
        "${phenotype}" \
        "${new_min_maf}"
    elif [ "${cluster}" = "sge" ]; then
     qsub -N "${slurm_jname}" \
        -o "${slurm_lname}.log" \
        -e "${slurm_lname}.errors.log" \
        -P lindgren.prjc \
        -wd $(pwd) \
        -t ${tasks}  \
        -q "${sge_queue}" \
        -pe shmem ${slurm_shmem} \
        "${bash_script}" \
        "${in_gmat}" \
        "${in_var}" \
        "${intervals}" \
        "${out_prefix}" \
        "${P_cutoff}" \
        "${max_iter}" \
        "${min_mac}" \
        "${grm_mtx}" \
        "${grm_sam}" \
        "${phenotype}" \
        "${new_min_maf}"
      fi
  else
    >&2 echo "${intervals} (interval) does not exist. Exiting.."
  fi

}

#submit_cts_analysis "pLoF_damaging_missense"
submit_binary_analysis "pLoF_damaging_missense"




