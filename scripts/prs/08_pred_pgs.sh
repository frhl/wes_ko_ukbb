#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=pred_pgs
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/pred_pgs.log
#SBATCH --error=logs/pred_pgs.errors.log
#SBATCH --open-mode=append
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-100
#SBATCH# --begin=now+8hour
#
#$ -N pred_pgs
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/pred_pgs.log
#$ -e logs/pred_pgs.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qe
#$ -t 1-100
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly bash_script="scripts/prs/_pred_pgs.sh"
readonly rscript="scripts/prs/08_pred_pgs.R"
readonly clean_script="scripts/prs/_prs_clean.sh"
readonly aggr_script="scripts/prs/_prs_aggr.sh"

readonly pred_dir="data/prs/hapmap/ukb_500k_hm3/validation"
readonly ld_dir="data/prs/hapmap/ld/matrix_unrel_kin"
readonly ldsc_dir="data/prs/ldsc"
readonly pheno_dir="data/phenotypes"
readonly beta_dir="data/prs/weights"
readonly out_dir="data/prs/scores_new/by_chrom"
readonly mrg_dir="data/prs/scores_new"

readonly cluster=$( get_current_cluster )
readonly index=$( get_array_task_id )

readonly file_binary="${pheno_dir}/dec22_phenotypes_binary_500k.tsv.gz"
readonly pheno_list_binary="${pheno_dir}/dec22_phenotypes_binary_200k_header.tsv"
readonly phenotype_binary=$( sed "${index}q;d" ${pheno_list_binary} )

# ldpred2 parameters
readonly impute="mean2"

mkdir -p ${out_dir}
mkdir -p ${mrg_dir}

predict_prs()
{ 
  local nslots=${1}
  local phenotype=${2}
  local pred="${pred_dir}/ukb_hapmap_500k_eur_chrCHR.bed"
  local ldsc="${ldsc_dir}/ldsc_${phenotype}.rds"
  local betas="${beta_dir}/weights.qc.${phenotype}.txt.gz"
  local out_prefix="${out_dir}/pgs.${phenotype}.chrCHR"
  local merged="${mrg_dir}/${phenotype}_pgs_chrom.txt.gz"

  local qsub_fit="_pgs_${phenotype}"
  local qsub_aggr="_aggr_${phenotype}"
  local qsub_clean="_clean_${phenotype}"

  # submit phenotype
  if [ ! -z ${phenotype} ]; then
    if [ -f "${betas}" ]; then
      if [ ! -f "${merged}" ]; then
          set_up_rpy
          fit_pgs
          aggr_pgs
          #clean_pgs
      else
        >&2 echo "${merged} already exists. Skipping.."
      fi
    fi
  fi
}


fit_pgs()
{
  readonly prs_jname="${qsub_fit}"
  readonly prs_lname="_pred_pgs"
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
      "${betas}" \
      "${impute}" \
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
      "${betas}" \
      "${impute}" \
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
      --output="logs/${aggr_lname}.log" \
      --error="logs/${aggr_lname}.errors.log" \
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

predict_prs "2" "${phenotype_binary}"

