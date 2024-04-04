#!/usr/bin/env bash
#
# condition on additive contribution (012)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=cond_additive_recessive
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/cond_additive_recessive.log
#SBATCH --error=logs/cond_additive_recessive.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=11-331
#SBATCH --open-mode=append

set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly cluster=$( get_current_cluster)
readonly index=$( get_array_task_id )

readonly rscript="scripts/saige/_check_prs_ok.R"

readonly curwd=$(pwd)
readonly vcf_dir="data/conditional/dominance/combine_encodings"
readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"

readonly grm_dir="data/saige/grm/input/dnanexus"
readonly grm_mtx="${grm_dir}/ukb_eur_200k_grm_fitted_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"
readonly plink_file="${grm_dir}/ukb_eur_200k_grm_grch38_rv_merged"

readonly spa_script="scripts/conditional/dominance/_cond_additive_recessive.sh"
readonly merge_script="scripts/_spa_merge.sh"
readonly in_prefix="ukb_eur_wes_200k"

# list of genes that passes significance cutoffs
readonly sig_genes_dir="data/conditional/combined/sig_genes"
readonly sig_genes="${sig_genes_dir}/sig_genes_after_sig_prs_176k.txt.gz"

# list of additive markers to condition on
readonly cond_additive_dir="data/conditional/dominance/combine_encodings"
readonly cond_additive_file="${cond_additive_dir}/ukb_eur_wes_200k_chrCHR_pLoF_damaging_missense.additive.txt"

# list of common markers to condition on
readonly cond_common_dir="data/conditional/common/combined"
readonly cond_common_file="${cond_common_dir}/ukb_eur_wes_200k_chrCHR_pLoF_damaging_missense_w_common_markers.txt"

# what categories should be included downstream?
readonly cond_cat="additive,common"

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

    local step1_dir="data/saige/output/${trait}/step1"
    local step2_dir="data/saige/output/${trait}/step2_dominance/min_mac${min_mac}"
    local in_vcf="${vcf_dir}/${in_prefix}_chrCHR_${annotation}.vcf.gz"
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
        >&2 echo "PRS enabled"
      else
        >&2 echo "Using without PRS."
      fi
    fi

    # setup paths to variants in genes by phenotypes
    #local markers_rare_by_gene="${markers_rare_by_gene_dir}/${in_prefix}_${phenotype}_${annotation}.txt.gz"

    if [ ! -f "${out_mrg}" ]; then
      local qsub_spa_name="spa_${phenotype}_${annotation}"
      local qsub_merge_name="_mrg_${phenotype}_${annotation}" 
      submit_spa_job
      #submit_merge_job
    else
      >&2 echo "Phenotype ${phenotype} with annotation ${annotation} already exists! Skipping.."
    fi
  else
    >&2 echo "No phenotype at index ${index}. Exiting.."
  fi
}



submit_spa_job() {
  mkdir -p ${step2_dir}
  local slurm_tasks="${tasks}"
  local slurm_jname="_${phenotype}_${annotation}"
  local slurm_lname="logs/_cond_additive_recessive"
  local slurm_project="${project}"
  local slurm_queue="short"
  local slurm_nslots="1"
  readonly spa_jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname}.log" \
    --error="${slurm_lname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --open-mode="append" \
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
    "${sig_genes}" \
    "${cond_additive_file}" \
    "${cond_common_file}" \
    "${cond_cat}" )
  echo "Submitted ${slurm_jname} ${spa_jid} "
}


submit_merge_job()
{
  local remove_by_chr="Y"
  local slurm_jname="_m_${phenotype}_${annotation}"
  local slurm_lname="logs/_mrg_cond_additive_recessive"
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
readonly use_prs="1"
readonly min_mac=4
readonly tasks=1-22 # 1-22
readonly project="lindgren.prj"

# cts traits
#submit_spa_cts_with_csqs "pLoF_damaging_missense"
submit_spa_binary_with_csqs "pLoF_damaging_missense"





