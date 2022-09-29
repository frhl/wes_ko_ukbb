#!/usr/bin/env bash
#
# Condition on common markers, all rare markers and PRS, to see if any
# associations are still significant. These would be driven by Compound heterozygotes.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=brute_force_cond
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/brute_force_cond.log
#SBATCH --error=logs/brute_force_cond.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=49

set -o errexit
set -o nounset


module purge
source utils/bash_utils.sh

readonly curwd=$(pwd)
readonly vcf_dir="data/conditional/rare/combined"
readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"

readonly grm_dir="data/saige/grm/input"
readonly grm_mtx="${grm_dir}/211102_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"

readonly spa_script="scripts/conditional/combined/_brute_force_cond.sh"
readonly merge_script="scripts/_spa_merge.sh"
readonly in_prefix="ukb_eur_wes_200k"

readonly markers_rare_by_gene_dir="data/conditional/rare/combined/genes/min_mac4"
readonly cond_rare_dir="data/conditional/rare/combined"

readonly cond_rare_file="${cond_rare_dir}/ukb_eur_wes_200k_chrCHR_maf0to5e-2_pLoF_damaging_missense_markers.txt.gz"
readonly markers_rare_ac="${cond_rare_dir}/ukb_eur_wes_200k_chrCHR_maf0to5e-2_pLoF_damaging_missense_AC.txt.gz"
readonly markers_rare_hash="${cond_rare_dir}/ukb_eur_wes_200k_chrCHR_maf0to5e-2_pLoF_damaging_missense_hash.txt.gz"

readonly cond_common_dir="data/conditional/common/combined/final"
readonly cond_common_file="${cond_common_dir}/ukb_eur_wes_200k_chrCHR_maf0to5e-2_pLoF_damaging_missense_w_common_markers.txt"

# what categories should be included downstream?
readonly cond_cat="pLoF,damaging_missense,common"

readonly index="${SLURM_ARRAY_TASK_ID}"

submit_spa_binary_with_csqs()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
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

    local step1_dir="data/saige/output/${trait}/step1"
    local step2_dir="data/saige/output/${trait}/step2_brute_force/min_mac${min_mac}"
    local in_vcf="${vcf_dir}/${in_prefix}_chrCHR_${maf}_${annotation}.vcf.bgz"
    mkdir -p ${step2_dir}

    local in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
    local in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
    local out_prefix="${step2_dir}/${in_prefix}_chrCHR_${maf}_${phenotype}_${annotation}"
    local out_mrg="${step2_dir}/${in_prefix}_${maf}_${phenotype}_${annotation}.txt.gz"

   if [ "${use_prs}" -eq "1" ]; then
      local in_gmat_prs="${step1_dir}/ukb_wes_200k_${phenotype}_chrCHR.rda"
      local in_var_prs="${step1_dir}/ukb_wes_200k_${phenotype}_chrCHR.varianceRatio.txt"
      if [ -f "${in_gmat_prs/CHR/21}" ] & [ -f "${in_var_prs/CHR/21}" ]; then
        local in_gmat=${in_gmat_prs}
        local in_var=${in_var_prs}
        local out_prefix="${step2_dir}/${in_prefix}_chrCHR_${maf}_${phenotype}_${annotation}_locoprs"
        local out_mrg="${step2_dir}/${in_prefix}_${maf}_${phenotype}_${annotation}_locoprs.txt.gz"
      else
        >&2 echo "Saige NULL (PRS) ${in_gmat_prs}/${in_var_prs} does not exist. Using without PRS."
      fi
    fi

    # setup paths to variants in genes by phenotypes
    local markers_rare_by_gene="${markers_rare_by_gene_dir}/${in_prefix}_${maf}_${phenotype}_${annotation}.txt.gz"

    if [ ! -f "${out_mrg}" ]; then
      submit_spa_job
      submit_merge_job
    else
      >&2 echo "Phenotype ${phenotype} with annotation ${annotation} already exists! Skipping.."
    fi
  else
    >&2 echo "No phenotype at index ${SLURM_ARRAY_TASK_ID}. Exiting.."
  fi
}


submit_spa_job() {
  mkdir -p ${step2_dir}
  local slurm_tasks="${tasks}"
  local slurm_jname="_${phenotype}_${annotation}"
  local slurm_lname="${out_prefix}"
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
    "${markers_rare_ac}" \
    "${markers_rare_hash}" \
    "${markers_rare_by_gene}" \
    "${markers_rare_cond_min_mac}" \
    "${cond_rare_file}" \
    "${cond_common_file}" \
    "${cond_cat}" )
  echo "Submitted brute force spa (jid=${spa_jid})"
}


submit_merge_job()
{
  local remove_by_chr="Y"
  local slurm_jname="_m_${phenotype}_${annotation}"
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
  echo "Submitted brute force merge (jid=${merge_jid})"
}


# parameters
readonly markers_rare_cond_min_mac=4
readonly use_prs="1"
readonly min_mac=4
readonly tasks=1-22
readonly project="lindgren.prj"
readonly queue="short"
readonly nslots=2

# Binary traits
maf="maf0to5e-2"

# cts traits
#submit_spa_cts_with_csqs "pLoF_damaging_missense"
submit_spa_binary_with_csqs "pLoF_damaging_missense"





