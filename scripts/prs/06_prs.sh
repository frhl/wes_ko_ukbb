#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=prs
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/prs.log
#SBATCH --error=logs/prs.errors.log
#SBATCH --open-mode=append
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1
#SBATCH --array=6
#
#
#$ -N prs
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prs.log
#$ -e logs/prs.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1-320 
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly bash_script="scripts/prs/_prs.sh"
readonly rscript="scripts/prs/06_prs.R"
readonly rscript_ldsc="scripts/prs/_check_ldsc.R"
readonly clean_script="scripts/prs/_prs_clean.sh"
readonly aggr_script="scripts/prs/_prs_aggr.sh"

readonly ldsc_dir="data/prs/ldsc"
readonly pred_dir="data/prs/hapmap/ukb_500k/validation"
readonly ld_dir="data/prs/hapmap/ld/matrix_unrel_kin"
readonly pheno_dir="data/phenotypes"
readonly out_dir="data/prs/scores/auto"
readonly mrg_dir="data/prs/scores"

# do not run files that have h2 estimates
# above the given p-value cutoff (nominal).
readonly ldsc_pvalue_cutoff="0.05"

readonly cluster=$( get_current_cluster )
readonly index=$( get_array_task_id )

readonly file_cts="${pheno_dir}/curated_covar_phenotypes_cts.tsv.gz"
readonly pheno_list_cts="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
readonly phenotype_cts=$( sed "${index}q;d" ${pheno_list_cts} )

readonly file_binary="${pheno_dir}/dec22_phenotypes_binary_500k.tsv.gz"
readonly pheno_list_binary="${pheno_dir}/dec22_phenotypes_binary_200k_header.tsv"
readonly phenotype_binary=$( sed "${index}q;d" ${pheno_list_binary} )

# ldpred2 parameters
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
  local merged="${mrg_dir}/${phenotype}_pgs.txt.gz"

  local qsub_fit="_prs_${phenotype}"
  local qsub_aggr="_aggr_${phenotype}"
  local qsub_clean="_clean_${phenotype}"

  # submit phenotype
  if [ ! -z ${phenotype} ]; then
    if [ ! -f "${merged}" ]; then
      # check that p-value passes thresholds
      set_up_rpy
      local pass_qc=$( Rscript ${rscript_ldsc} --ldsc ${ldsc} --ldsc_pvalue_cutoff ${ldsc_pvalue_cutoff} )
      if [ "${pass_qc}" = "1" ]; then
        fit_pgs
        aggr_pgs
        clean_pgs
      else
        >&2 echo "${phenotype} does not pass LDSC QC."
      fi
    else
      >&2 echo "${merged} already exists. Skipping.."
    fi
    # remove disk backing files
  fi

}


fit_pgs()
{
  readonly prs_jname="${qsub_fit}"
  readonly prs_lname="_prs"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local sge_queue="short.qc"
  local slurm_tasks="${tasks}"
  local slurm_nslots="${nslots}"
  fit_pgs_jid=""
  if [ "${cluster}" = "slurm" ]; then
    fit_pgs_jid=$( sbatch \
      --account="${slurm_project}" \
      --job-name="${prs_jname}" \
      --output="logs/${prs_lname}.log" \
      --error="logs/${prs_lname}.errors.log" \
      --chdir="$(pwd)" \
      --partition="${slurm_queue}" \
      --cpus-per-task="${slurm_nslots}" \
      --open-mode=append \
      --array=${slurm_tasks} \
      --parsable \
      "${bash_script}" \
      "${rscript}" \
      "${pred}" \
      "${ldsc}" \
      "${ld_dir}" \
      "${method}" \
      "${impute}" \
      "${ldsc_pvalue_cutoff}" \
      "${out_prefix}" )
  elif [ "${cluster}" = "sge" ]; then
    qsub -N "${prs_jname}" \
      -t ${slurm_tasks} \
      -q ${sge_queue} \
      -o "logs/${prs_lname}.log" \
      -e "logs/${prs_lname}.errors.log" \
      -P lindgren.prjc \
      -wd $(pwd) \
      -pe shmem "${slurm_nslots}" \
      "${bash_script}" \
      "${rscript}" \
      "${pred}" \
      "${ldsc}" \
      "${ld_dir}" \
      "${method}" \
      "${impute}" \
      "${ldsc_pvalue_cutoff}" \
      "${out_prefix}"
  else
    >&2 echo "${cluster} is not valid!"
  fi

}


aggr_pgs()
{
  readonly aggr_jname="${qsub_aggr}"
  readonly aggr_lname="_aggr_prs"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local sge_queue="short.qc"
  local slurm_tasks="1"
  local slurm_nslots="1"
  if [ "${cluster}" = "slurm" ]; then
    readonly aggr_pgs_jid=$( sbatch \
      --account="${slurm_project}" \
      --job-name="${aggr_jname}" \
      --output="logs/${aggr_jname}.log" \
      --error="logs/${aggr_jname}.errors.log" \
      --chdir="$(pwd)" \
      --partition="${slurm_queue}" \
      --cpus-per-task="${slurm_nslots}" \
      --array=${slurm_tasks} \
      --open-mode=append \
      --dependency="afterok:${fit_pgs_jid}" \
      --parsable \
      "${aggr_script}" \
      "${phenotype}" \
      "${out_dir}" \
      "${mrg_dir}" )
  elif [ "${cluster}" = "sge" ]; then
    qsub -N "${qsub_aggr}" \
      -q test.qc \
      -P lindgren.prjc \
      -o "logs/${aggr_lname}.log" \
      -e "logs/${aggr_lname}.errors.log" \
      -pe shmem 1 \
      -hold_jid "${qsub_fit}" \
      -wd $(pwd) \
      "${aggr_script}" \
      "${phenotype}" \
      "${out_dir}" \
      "${mrg_dir}"
  else
    >&2 echo "${cluster} is not valid"
  fi

}


clean_pgs()
{
  readonly clean_jname="${qsub_clean}"
  readonly clean_lname="_prs_clean"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local slurm_tasks="${tasks}"
  local slurm_nslots="1"
  if [ "${cluster}" = "slurm" ]; then
    readonly clean_pgs_jid=$( sbatch \
      --account="${slurm_project}" \
      --job-name="${clean_jname}" \
      --output="logs/${clean_lname}.log" \
      --error="logs/${clean_lname}.errors.log" \
      --chdir="$(pwd)" \
      --partition="${slurm_queue}" \
      --cpus-per-task="${slurm_nslots}" \
      --array=${slurm_tasks} \
      --dependency="aftercorr:${fit_pgs_jid}" \
      --open-mode=append \
      --parsable \
      "${clean_script}" \
      "${pred}" \
      "${out_prefix}" )
  elif [ "${cluster}" = "sge" ]; then
    qsub -N "${qsub_clean}" \
      -t ${tasks} \
      -q test.qc \
      -P lindgren.prjc \
      -o "logs/${clean_lname}.log" \
      -e "logs/${clean_lname}.errors.log" \
      -wd $(pwd) \
      -pe shmem 1 \
      -hold_jid_ad "${qsub_fit}" \
      "${clean_script}" \
      "${pred}" \
      "${out_prefix}" 
  else
    >&2 echo "${cluster} is not valid"
  fi

}

# parameters
readonly queue="short"
readonly project="lindgren.prj"
readonly tasks=1-22

submit_ldpred2 "auto" "2" "${phenotype_binary}"
#submit_ldpred2 "auto" "6" "${phenotype_cts}_int"
#submit_ldpred2 "auto" "6" "${phenotype_cts}"

