#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=spa_test_group
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/spa_test_group.log
#SBATCH --error=logs/spa_test_group.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-22
#SBATCh --dependency="afterok:39656974"

set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh
source utils/qsub_utils.sh

#readonly vcf_dir="data/mt/prefilter/pp90" # with standard PP cutoff
readonly vcf_dir="data/mt/prefilter/no_pp_cutoff/old" # without PP cutoff
readonly pheno_dir="data/phenotypes"

readonly plink_dir="data/saige/grm/input"
readonly grm_dir="data/saige/grm/input/dnanexus"

readonly grm_mtx="${grm_dir}/ukb_eur_200k_grm_fitted_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"
readonly plink_file="${grm_dir}/ukb_eur_200k_grm_grch38_rv_merged"

readonly vep_version="95"
readonly saige_group_file_dir="data/vep/vep${vep_version}/saige_group"
readonly saige_group_file="${saige_group_file_dir}/UKB.chrCHR.exome_array.variants_only.vep${vep_version}.csqs.worst_csq_by_gene_canonical.original.saige.txt.gz"

readonly rscript="scripts/saige/_check_prs_ok.R"
readonly spa_script="scripts/saige/_spa_test_group.sh"
readonly merge_script="scripts/saige/_spa_merge_group.sh"
readonly in_prefix="ukb_eur_wes_200k"

readonly cluster=$( get_current_cluster)
readonly index=$( get_array_task_id )
readonly curwd=$(pwd)

submit_spa_binary_with_csqs()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  #local pheno_list="${pheno_dir}/dec22_phenotypes_binary_200k_header.tsv"
  local pheno_list="${pheno_dir}/dec22_phenotypes_binary_200k_header_additive_to_try.txt"
  local phenotype=$( sed "${index}q;d" ${pheno_list} )
  submit_spa_with_csqs "${annotation}" "${phenotype}" "binary"
}

submit_spa_with_csqs()
{
  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local trait=${3?Error: Missing arg3 (trait)}
  local annotation_label="$(echo $annotation | sed -e 's/,/_/g')"
  if [ ! -z ${phenotype} ]; then

    local step1_dir="data/saige/output/${trait}/step1"
    #local step2_dir="data/saige/output/${trait}/step2_saige_group/min_mac${min_mac}"
    local step2_dir="data/saige/output/${trait}/step2_saige_group_no_pp/min_mac${min_mac}"
    local in_vcf="${vcf_dir}/ukb_wes_union_calls_200k_chrCHR.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.vcf.bgz"
    mkdir -p ${step2_dir}

    local in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
    local in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
    local out_prefix="${step2_dir}/${in_prefix}_chrCHR_${phenotype}_pLoF_damaging_missense.vep${vep_version}"
    local out_mrg="${step2_dir}/${in_prefix}_${phenotype}_${annotation_label}.vep${vep_version}.txt.gz"

    if [ "${use_prs}" -eq "1" ]; then
      set_up_rpy
      local in_gmat_prs="${step1_dir}/ukb_wes_200k_${phenotype}_chrCHR.rda"
      local in_var_prs="${step1_dir}/ukb_wes_200k_${phenotype}_chrCHR.varianceRatio.txt"
      local prs_ok=$(Rscript ${rscript} --phenotype ${phenotype} --include_nominal_significant )
      if [ -f "${in_gmat_prs/CHR/21}" ] & [ -f "${in_var_prs/CHR/21}" ] & [ "${prs_ok}" -eq "1" ]; then
        local in_gmat=${in_gmat_prs}
        local in_var=${in_var_prs}
        local out_prefix="${step2_dir}/${in_prefix}_chrCHR_${phenotype}_${annotation_label}_locoprs.vep${vep_version}"
        local out_mrg="${step2_dir}/${in_prefix}_${phenotype}_${annotation_label}_locoprs.vep${vep_version}.txt.gz"
      else
        >&2 echo "Using without PRS."
      fi 
    fi

    if [ ! -f "${out_mrg}" ]; then
      local slurm_spa_name="spa_${phenotype}_${annotation_label}"
      local slurm_merge_name="_mrg_${phenotype}_${annotation_label}"
      jid=$(submit_spa_job)
      submit_merge_job ${jid} ## need to modify merge
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
  local slurm_jname="${slurm_spa_name}"
  local slurm_lname="logs/_spa_test_group"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local slurm_nslots="${nslots}"
  local spa_jid=$( sbatch \
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
    "${saige_group_file}" \
    "${annotation}" \
    "${max_maf_in_group_test}" \
    "${out_prefix}" )
  echo ${spa_jid}
}


submit_merge_job()
{
  local waitfor="${1}"
  local remove_by_chr="Y"
  local slurm_jname="${slurm_merge_name}"
  local slurm_lname="logs/_spa_merge"
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
    --dependency="afterok:${waitfor}" \
    --open-mode="append" \
    "${merge_script}" \
    "${out_prefix}" \
    "${out_mrg}" \
    "${remove_by_chr}")
}

# parameters
readonly use_prs="0"
readonly min_mac=4
readonly project="lindgren.prj"
readonly tasks="1-22"
readonly queue="short"
readonly nslots=2
readonly max_maf_in_group_test="0.05"

# pLoF_damaging_missense (combinesi pLoF and damaging_missense into a single category)  
# pLoF,damaging_missense (runs pLoF and damaging_missense seperately)

submit_spa_binary_with_csqs "pLoF_damaging_missense"





