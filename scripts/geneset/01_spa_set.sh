#!/usr/bin/env bash
#
# @description run SAIGE-GENE+ using current variant annotations.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=spa_set
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/spa_set.log
#SBATCH --error=logs/spa_set.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-22

# 12,18

set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh

#readonly vcf_dir="data/knockouts/alt"
readonly vcf_dir="data/mt/annotated"
readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"
readonly grm_dir="data/saige/grm/input"

readonly spa_script="scripts/geneset/_spa_set.sh"
readonly merge_script="scripts/_spa_merge.sh"
readonly in_prefix="ukb_eur_wes_union_calls_200k"
readonly in_vcf="${vcf_dir}/${in_prefix}_chrCHR.vcf.bgz"

readonly grm_mtx="${grm_dir}/211102_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"

readonly group_dir="data/mt/vep/worst_csq_by_gene_canonical"
readonly group="${group_dir}/ukb_eur_wes_union_calls_200k_chrCHR.saige"
readonly index="${SLURM_ARRAY_TASK_ID}"

submit_spa_set_binary()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
  local phenotype=$( sed "${index}q;d" ${pheno_list} )
  submit_spa_pair "${annotation}" "${phenotype}" "binary"
}

submit_spa_set_cts()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( sed "${index}q;d" ${pheno_list} )
  submit_spa_pair "${annotation}" "${phenotype}" "cts"
}

submit_spa_pair()
{

  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local trait=${3?Error: Missing arg3 (trait)}

  local step1_dir="data/saige/output/${trait}/step1"
  local step2_dir="data/saige/output/${trait}/step2_set_AUG/min_mac${min_mac}"
  mkdir -p ${step2_dir}

  local in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
  local in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
  local out_prefix="${step2_dir}/${in_prefix}_chrCHR_${maf}_${phenotype}_${annotation}"
  local out_mrg="${step2_dir}/${in_prefix}_${maf}_${phenotype}_${annotation}.txt.gz"

  # setup PRS
  if [ "${use_prs}" -eq "1" ]; then
    local in_gmat_prs="${step1_dir}/ukb_wes_200k_${phenotype}_chrCHR.rda"
    local in_var_prs="${step1_dir}/ukb_wes_200k_${phenotype}_chrCHR.varianceRatio.txt"
    if [ -f "${in_gmat_prs/CHR/21}" ] & [ -f "${in_var_prs/CHR/21}" ]; then
      local in_gmat=${in_gmat_prs}
      local in_var=${in_var_prs}
      local out_prefix="${step2_dir}/${in_prefix}_chrCHR_${maf}_${phenotype}_${annotation}_locoprs"
      local out_mrg="${step2_dir}/${in_prefix}_${maf}_${phenotype}_${annotation}_locoprs.txt.gz"
    else
      >&2 echo "Saige NULL (PRS) ${in_gmat_prs} or ${in_var_prs} does not exist. Using without PRS."
    fi
  fi 

  # Submit scripts
  if [ -f "${in_gmat/CHR/21}" ] && [ -f "${in_var/CHR/21}" ]; then
    if [ ! -f "${out_mrg}" ]; then
      local slurm_spa_name="sspa_${phenotype}_${annotation}"
      local slurm_merge_name="_smrg_${phenotype}_${annotation}"  
      submit_spa_set_job
      submit_merge_job
    else
      >&2 echo "${out_mrg} already exists. Skipping.."
    fi
  else
    >&2 echo "${in_gmat} and/or ${in_var} does not exist!"
  fi

}


submit_spa_set_job() 
{
  mkdir -p ${step2_dir}
  local slurm_tasks="${tasks}"
  local slurm_jname="${slurm_spa_name}"
  local slurm_lname="logs/_spa_set_test"
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
    "${group}" \
    "${out_prefix}" )
  echo "Submitting ${phenotype} (JID=${spa_jid}).."
}


submit_merge_job()
{
  local remove_by_chr="Y"
  local slurm_jname="${slurm_merge_name}"
  local slurm_lname="logs/_spa_set_merge"
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

readonly maf="maf0to5e-2"
readonly project="lindgren.prj"
readonly min_mac=4
readonly tasks=1-2
readonly queue="short"
readonly nslots=1
readonly use_prs="1"

# cts traits
submit_spa_set_binary "pLoF_damaging_missense"
#submit_spa_set_cts "pLoF_damaging_missense"


