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
#SBATCH --array=1-5
#SBATCH --requeue


set -o errexit
set -o nounset

readonly curwd=$(pwd)
readonly bash_script="scripts/conditional/common/_spa_iter_common.sh"


# parameters
readonly min_maf=0.01
readonly min_mac=4
readonly max_iter=30
readonly P_cutoff="5e-6"

# directories and paths
readonly pheno_dir="data/phenotypes"
readonly interval_dir="data/conditional/common/intervals/min_mac${min_mac}"
readonly out_dir="data/conditional/common/spa_iter/"

readonly grm_dir="data/saige/grm/input/dnanexus"
readonly grm_mtx="${grm_dir}/ukb_eur_200k_grm_fitted_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"

readonly in_prefix="ukb_eur_wes_200k"

submit_binary_analysis()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/spiros_brava_phenotypes_binary_200k_header.tsv"
  local phenotype=$( sed "${SLURM_ARRAY_TASK_ID}q;d" ${pheno_list} )
  submit_cond_spa "${annotation}" "${phenotype}" "binary"
}

submit_cts_analysis()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( sed "${SLURM_ARRAY_TASK_ID}q;d" ${pheno_list} )
  submit_cond_spa "${annotation}" "${phenotype}" "cts"
}

submit_cond_spa()
{
  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local trait=${3?Error: Missing arg3 (trait)}

  local step1_dir="data/saige/output/${trait}/step1"
  local step2_dir="data/saige/output/${trait}/step2/minmac${min_mac}"
  local in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
  local in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
  local interval_vcf="${interval_dir}/${in_prefix}_${phenotype}_${annotation}.vcf.bgz"
  local out_prefix="${out_dir}/${in_prefix}_${phenotype}_${annotation}_cond"

  echo "interval_vcf: $( ls -l ${interval_vcf})"

  mkdir -p ${out_dir}
  if [ -f "${interval_vcf}" ]; then 
    readonly slurm_tasks="${SLURM_ARRAY_TASK_ID}"
    readonly slurm_jname="_cond_${phenotype}"
    readonly slurm_lname="${out_prefix}"
    readonly slurm_project="lindgren.prj"
    readonly slurm_queue="long"
    readonly slurm_shmem="1"
    set -x
    sbatch \
      --account="${slurm_project}" \
      --job-name="${slurm_jname}" \
      --output="${slurm_lname}.log" \
      --error="${slurm_lname}.errors.log" \
      --chdir="${curwd}" \
      --partition="${slurm_queue}" \
      --cpus-per-task="${slurm_shmem}" \
      --array=${slurm_tasks} \
      --parsable \
      "${bash_script}" \
      "${in_gmat}" \
      "${in_var}" \
      "${interval_vcf}" \
      "${out_prefix}" \
      "${P_cutoff}" \
      "${max_iter}" \
      "${min_mac}" \
      "${grm_mtx}" \
      "${grm_sam}" \
      "${phenotype}" \
      "${min_maf}"
    set +x
  else
    >&2 echo "${interval_vcf} (interval) does not exist. Exiting.."
  fi

}

#submit_cts_analysis "pLoF_damaging_missense"
submit_binary_analysis "pLoF_damaging_missense"




