#!/usr/bin/env bash
#
# @description: amalgamate variants by phase to infer knockouts by genes.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=extract_knockouts
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/extract_knockouts.log
#SBATCH --error=logs/extract_knockouts.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1
#SBATCH --array=21
#
#
#$ -N extract_knockouts
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/extract_knockouts.log
#$ -e logs/extract_knockouts.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1,5
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly curwd=$(pwd)
readonly spark_dir="data/tmp/spark"
readonly bash_script="scripts/_extract_knockouts.sh"
readonly hail_script="scripts/_write_gene_intervals.py"

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_dir="data/mt/prefilter/final_99"
readonly out_dir="data/knockouts/extracted/chr${chr}"
readonly in_prefix="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp99.maf0_005.mt"
readonly in_type="mt"

readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}"
readonly out_interval="${out_prefix}_interval.txt"

readonly queue="epyc"
readonly project="lindgren.prj"
readonly only_vcf=""

mkdir -p ${out_dir}

# create a file containing genes used in analysis
extract_genes() {
  if [ ! -f "${out_interval}" ]; then
    set_up_hail
    set_up_pythonpath_legacy
    python3 "${hail_script}" \
      --input_path ${in_prefix} \
      --input_type ${in_type} \
      --out_prefix ${out_interval}
 else
  >&2 echo "${out_interval} already exists. Skipping."
 fi
}


# Get number of genes to run
extract_genes
readonly n_genes=$( cat ${out_interval} | wc -l )
readonly array_id="1" #-${n_genes}"


submit_knockout_job() 
{
  # I/O
  local annotation=${1}
  local nslots=${2}
  local out_prefix_csqs="${out_prefix}_${annotation/,/_}"
  
  # slurm specific paramters 
  local slurm_jname="_c${chr}_extr_${annotation}"
  local slurm_lname="logs/_extract_knockouts"
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
      --array=${array_id} \
      --parsable \
      "${bash_script}" \
      "${in_prefix}" \
      "${in_type}" \
      "${annotation}" \
      "${out_interval}" \
      "${out_prefix_csqs}" \
      "${chr}"
  elif [ "${cluster}" = "sge" ]; then
    qsub -N "${slurm_jname}" \
      -o "${slurm_lname}.log" \
      -e "${slurm_lname}.errors.log" \
      -P lindgren.prjc \
      -wd $(pwd) \
      -t ${array_id} \
      -q "${sge_queue}" \
      -pe shmem ${slurm_nslots} \
      "${bash_script}" \
      "${in_prefix}" \
      "${in_type}" \
      "${annotation}" \
      "${out_interval}" \
      "${out_prefix_csqs}" \
      "${chr}"
  else
    >&2 echo "${cluster} is not a valid cluster."
  fi
}

submit_knockout_job "pLoF,damaging_missense" "1"



