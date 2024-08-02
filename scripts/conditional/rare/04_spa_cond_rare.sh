#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=spa_cond_rare
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/spa_cond_rare.log
#SBATCH --error=logs/spa_cond_rare.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1
#SBATCH --array=11-50

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly cluster=$( get_current_cluster)
readonly index=$( get_array_task_id )

readonly vcf_dir="data/conditional/rare/combined/mt"
readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"

readonly grm_dir="data/saige/grm/input/dnanexus"
readonly grm_mtx="${grm_dir}/ukb_eur_200k_grm_fitted_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"
readonly plink_file="${grm_dir}/ukb_eur_200k_grm_grch38_rv_merged"

readonly spa_script="scripts/conditional/rare/_spa_cond_rare.sh"
readonly merge_script="scripts/saige/_spa_merge.sh"
readonly in_prefix="ukb_eur_wes_200k"

# directory to markers by genes that are significant in primary analysis
readonly markers_by_gene_dir="data/conditional/rare/combined/genes/2024/min_mac4"

# directory to conditioning markers
readonly cond_rare_dir="data/conditional/rare/combined/mt"
readonly cond_rare_file="${cond_rare_dir}/ukb_eur_wes_200k_chrCHR_pLoF_damaging_missense_markers.txt.gz"

# path to file with allele count by phenotype (need to avoid conditioning on monomorphic SNPs)
readonly markers_ac="${cond_rare_dir}/ukb_eur_wes_200k_chrCHR_pLoF_damaging_missense_AC.txt.gz"
readonly markers_hash="${cond_rare_dir}/ukb_eur_wes_200k_chrCHR_pLoF_damaging_missense_hash.txt.gz"
# category group for markers to condition on (csv)
readonly cond_cat="pLoF,damaging_missense" 


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

    # setup input and outputs
    local step1_dir="data/saige/output/${trait}/step1"
    local step2_dir="data/saige/output/${trait}/step2_rare_cond/min_mac${min_mac}"
    local in_vcf="${vcf_dir}/${in_prefix}_chrCHR_${annotation}.vcf.bgz"
    mkdir -p ${step2_dir}

    # setup paths to saige step 1
    local in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
    local in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
    local out_prefix="${step2_dir}/${in_prefix}_chrCHR_${phenotype}_${annotation}"
    local out_mrg="${step2_dir}/${in_prefix}_${phenotype}_${annotation}.txt.gz"

    # setup paths to saige step 1 PRS scores
   if [ "${use_prs}" -eq "1" ]; then
      local in_gmat_prs="${step1_dir}/ukb_wes_200k_${phenotype}_chrCHR.rda"
      local in_var_prs="${step1_dir}/ukb_wes_200k_${phenotype}_chrCHR.varianceRatio.txt"
      if [ -f "${in_gmat_prs/CHR/21}" ] & [ -f "${in_var_prs/CHR/21}" ]; then
        local in_gmat=${in_gmat_prs}
        local in_var=${in_var_prs}
        local out_prefix="${step2_dir}/${in_prefix}_chrCHR_${phenotype}_${annotation}_locoprs"
        local out_mrg="${step2_dir}/${in_prefix}_${phenotype}_${annotation}_locoprs.txt.gz"
      else
        >&2 echo "Saige NULL (PRS) ${in_gmat_prs}/${in_var_prs} does not exist. Using without PRS."
      fi
    fi

    # setup paths to variants in genes by phenotype (We don't care about PRS here, since
    # this step is just for selecting variants within genes).
    #updated_phenotype=$(echo ${phenotype} | sed -e s/"_primary_care"/""/)
    local markers_by_gene="${markers_by_gene_dir}/${in_prefix}_${phenotype}_${annotation}.txt.gz"

    if [ -f ${markers_by_gene} ]; then
      if [ ! -f "${out_mrg}" ]; then
        local slurm_spa_name="spa_${phenotype}_${annotation}"
        local slurm_merge_name="_mrg_${phenotype}_${annotation}"
        submit_spa_job
        submit_merge_job
      else
        >&2 echo "Phenotype ${phenotype} with annotation ${annotation} already exists! Skipping.." 
      fi
    else
      >&2 echo "Markers by gene file does not exists (${markers_by_gene})"
    fi 
  else
    >&2 echo "No phenotype at index ${index}. Exiting.." 
  fi 
}

submit_spa_job() {
  mkdir -p ${step2_dir}
  local slurm_tasks="${tasks}"
  local slurm_jname="${slurm_spa_name}"
  local slurm_lname="${out_prefix}"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local sge_queue="epyc.qc"
  local slurm_nslots="${nslots}"
  readonly spa_jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname}.log" \
    --error="${slurm_lname}.errors.log" \
    --chdir="$(pwd)" \
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
    "${markers_ac}" \
    "${markers_hash}" \
    "${markers_by_gene}" \
    "${markers_cond_min_mac}" \
    "${out_prefix}" \
    "${cond_rare_file}" \
    "${cond_cat}" )
}


submit_merge_job() {
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
    --chdir="$(pwd)" \
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
readonly markers_cond_min_mac=4  #3
readonly use_prs="0"
readonly min_mac=4
readonly tasks=1-22
readonly queue="epyc"
readonly project="lindgren.prj"
readonly nslots=4



# cts traits
#submit_spa_cts_with_csqs "pLoF_damaging_missense"
submit_spa_binary_with_csqs "pLoF_damaging_missense"



