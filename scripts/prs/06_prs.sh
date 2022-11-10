#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=prs
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/prs.log
#SBATCH --error=logs/prs.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=17-18
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly bash_script="scripts/prs/_prs.sh"
readonly rscript="scripts/prs/06_prs.R"
readonly clean_script="scripts/prs/_prs_clean.sh"
readonly aggr_script="scripts/prs/_prs_aggr.sh"

readonly ldsc_dir="data/prs/ldsc"
readonly pred_dir="data/prs/hapmap/ukb_500k/validation"
readonly ld_dir="data/prs/hapmap/ld/matrix_unrel_kin"
readonly pheno_dir="data/phenotypes"
readonly out_dir="data/prs/scores/auto"
readonly mrg_dir="data/prs/scores"

readonly index=${SLURM_ARRAY_TASK_ID}

readonly file_cts="${pheno_dir}/curated_covar_phenotypes_cts.tsv.gz"
readonly pheno_list_cts="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
readonly phenotype_cts=$( sed "${index}q;d" ${pheno_list_cts} )

readonly file_binary="${pheno_dir}/curated_covar_phenotypes_binary.tsv.gz"
readonly pheno_list_binary="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
readonly phenotype_binary=$( sed "${index}q;d" ${pheno_list_binary} )
readonly impute="mean2"

mkdir -p ${out_dir}

submit_ldpred2()
{ 
  local method=${1}
  local nslots=${2}
  local phenotype=${3}
  local pred="${pred_dir}/ukb_hapmap_500k_eur_chrCHR.bed"
  local ldsc="${ldsc_dir}/ldsc_${phenotype}.rds"
  local out_prefix="${out_dir}/prs_${method}_${phenotype}_chrCHR"

  local qsub_fit="_prs_${phenotype}"
  local qsub_aggr="_aggr_${phenotype}"
  local qsub_clean="_clean_${phenotype}"

  if [ ! -z ${phenotype} ]; then
    # fit actual pgs
    fit_pgs
    # aggregate into matrix
    aggr_pgs
    # remove disk backing files
    clean_pgs
  fi

}


fit_pgs()
{
  local slurm_jname="_pgs"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local slurm_tasks="${tasks}"
  local slurm_nslots="${nslots}"
  readonly fit_pgs_jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_jname}.log" \
    --error="${slurm_jname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --array=${slurm_tasks} \
    --parsable \
    "${bash_script}" \
    "${rscript}" \
    "${pred}" \
    "${ldsc}" \
    "${ld_dir}" \
    "${method}" \
    "${impute}" \
    "${out_prefix}" )
}


aggr_pgs()
{
  local slurm_jname="${qsub_aggr}"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local slurm_tasks="1"
  local slurm_nslots="1"
  readonly aggr_pgs_jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_jname}.log" \
    --error="${slurm_jname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --array=${slurm_tasks} \
    --dependency="afterok:${fit_pgs_jid}" \
    --parsable \
    "${aggr_script}" \
    "${phenotype}" \
    "${out_dir}" \
    "${mrg_dir}" )
}


clean_pgs()
{
  local slurm_jname="${qsub_clean}"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local slurm_tasks="${tasks}"
  local slurm_nslots="1"
  readonly clean_pgs_jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_jname}.log" \
    --error="${slurm_jname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --array=${slurm_tasks} \
    --dependency="aftercorr:${fit_pgs_jid}" \
    --parsable \
    -hold_jid_ad "_prs_${phenotype}" \
    "${clean_script}" \
    "${pred}" \
    "${out_prefix}" )
}

# parameters
readonly queue="short"
readonly project="lindgren.prj"
readonly tasks=1-22

submit_ldpred2 "auto" "6" "${phenotype_cts}_int"
submit_ldpred2 "auto" "6" "${phenotype_cts}"

