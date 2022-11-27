#!/usr/bin/env bash
#
# @description: amalgamate variants by phase to infer knockouts by genes.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=knockouts
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/knockouts.log
#SBATCH --error=logs/knockouts.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --requeue
#
#
#$ -N knockouts
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/knockouts.log
#$ -e logs/knockouts.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1
#$ -V
#$ -hold_jid 79672533

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly curwd=$(pwd)
readonly spark_dir="data/tmp/spark"
readonly bash_script="scripts/_encode_knockouts.sh"

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_dir="data/mt/prefilter/final_99"
readonly out_dir="data/knockouts/alt"
readonly in_prefix="${in_dir}/ukb_wes_union_calls_200k_chrCHR.loftee.worst_csq_by_gene_canonical.pp99.maf0_005.mt"
readonly in_type="mt"

readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chrCHR"
readonly out_type="vcf"

# Note: ~24 slots are needed for running chr1. 
# Note: long queue may be required for chr1.
readonly queue="short"
readonly project="lindgren.prj"

# should only VCF be produced?
readonly only_vcf=""

mkdir -p ${out_dir}


submit_knockout_job() 
{
  # I/O
  local annotation=${1}
  local nslots=${2}
  local aggr_method=${3}
  local out_prefix_csqs="${out_prefix}_${annotation/,/_}"
  local out_checkpoint="${out_prefix_csqs}_checkpoint.mt"
  
  # slurm specific paramters 
  local slurm_jname="_c${chr}_ko_${annotation}"
  local slurm_lname="logs/_knockouts"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local sge_queue="short.qe"
  local slurm_nslots="${nslots}"
  if [ "${cluster}" = "slurm" ]; then
    sbatch \
      --account="${slurm_project}" \
      --job-name="${slurm_jname}" \
      --output="${slurm_lname}.log" \
      --error="${slurm_lname}.errors.log" \
      --chdir="${curwd}" \
      --partition="${slurm_queue}" \
      --cpus-per-task="${slurm_nslots}" \
      --array=${task_id} \
      --parsable \
      "${bash_script}" \
      "${in_prefix}" \
      "${in_type}" \
      "${annotation}" \
      "${only_vcf}" \
      "${aggr_method}" \
      "${out_prefix_csqs}" \
      "${out_type}" 
  elif [ "${cluster}" = "sge" ]; then
    qsub -N "${slurm_jname}" \
      -o "${slurm_lname}.log" \
      -e "${slurm_lname}.errors.log" \
      -P lindgren.prjc \
      -wd $(pwd) \
      -t ${task_id} \
      -q "${sge_queue}" \
      -pe shmem ${slurm_nslots} \
      "${bash_script}" \
      "${in_prefix}" \
      "${in_type}" \
      "${annotation}" \
      "${only_vcf}" \
      "${aggr_method}" \
      "${out_prefix_csqs}" \
      "${out_type}" 
  else
    >&2 echo "${cluster} is not a valid cluster."
  fi
  # clean up after checkpints when
  if [ -f "${out_checkpoint}" ]; then
    rm -rf ${out_checkpoint}
  fi
}

submit_knockout_job "pLoF,damaging_missense" "30" "collect"
#submit_knockout_job "pLoF" "16" "collect"
#submit_knockout_job "damaging_missense" "30" "collect"
#submit_knockout_job "synonymous" "30" "collect"



