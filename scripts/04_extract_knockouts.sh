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
#SBATCH --array=4-5
#
#
#$ -N extract_knockouts
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/extract_knockouts.log
#$ -e logs/extract_knockouts.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1-3
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly curwd=$(pwd)
readonly spark_dir="data/tmp/spark"
readonly bash_script="scripts/_extract_knockouts.sh"
readonly merge_script="scripts/_merge_knockouts.sh"
readonly hail_script="scripts/_write_gene_intervals.py"

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_dir="data/mt/prefilter/final_90"
readonly merge_dir="data/knockouts/extracted/pp90"
readonly out_dir="data/knockouts/extracted/pp90/chr${chr}"
# in parameters
readonly in_prefix="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.mt"
readonly in_type="mt"
# prefix for indiviual genes and final merged file
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}"
readonly out_merge="${merge_dir}/ukb_eur_wes_200k_chr${chr}"
# interval file to keep track of genes to assess
readonly out_interval="${out_prefix}_interval.txt"

# slurm/sge paramters
readonly sge_project="lindgren.prjc"
readonly slurm_project="lindgren.prj"
readonly slurm_queue="epyc"
readonly sge_queue="short.qc" 

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
readonly array_id="1-${n_genes}"

submit_merge_job()
{
  echo "Submitting merge job."
  local regex_prefix="${1}"
  local outfile="${2}"
  local sge_dependency="${3}"
  local slurm_dependency="${4}"
  local jname="_c${chr}_mrg_ko"
  local lname="logs/_merge_knockouts"
  local nslots="1"
  if [ "${cluster}" = "slurm" ]; then
    sbatch \
      --account="${slurm_project}" \
      --job-name="${jname}" \
      --output="${lname}.log" \
      --error="${lname}.errors.log" \
      --chdir="${curwd}" \
      --partition="${slurm_queue}" \
      --cpus-per-task="${nslots}" \
      --dependency="${slurm_dependency}" \
      "${merge_script}" \
      "${regex_prefix}" \
      "${out_interval}" \
      "${outfile}" 
  elif [ "${cluster}" = "sge" ]; then
     qsub -N "${jname}" \
      -o "${lname}.log" \
      -e "${lname}.errors.log" \
      -P ${sge_project} \
      -wd $(pwd) \
      -q "${sge_queue}" \
      -pe shmem ${nslots} \
      -hold_jid ${sge_dependency} \
      "${merge_script}" \
      "${regex_prefix}" \
      "${out_interval}" \
      "${outfile}"
  fi
}


submit_knockout_job() 
{
  # I/O
  local annotation=${1}
  local nslots=${2}
  local out_prefix_csqs="${out_prefix}_${annotation/,/_}"
  local out_prefix_merge="${out_merge}_${annotation/,/_}"
  local merge_prefix_regex="${out_prefix_csqs}_ENSG"

  # slurm specific paramters 
  local ko_jname="_c${chr}_extr_${annotation}"
  local ko_lname="logs/_extract_knockouts"
  local slurm_nslots="${nslots}"
  local slurm_jid="na"
  if [ "${cluster}" = "slurm" ]; then
    local slurm_jid=$( sbatch \
      --account="${slurm_project}" \
      --job-name="${ko_jname}" \
      --output="${ko_lname}.log" \
      --error="${ko_lname}.errors.log" \
      --chdir="${curwd}" \
      --partition="${slurm_queue}" \
      --cpus-per-task="${nslots}" \
      --array=${array_id} \
      --parsable \
      "${bash_script}" \
      "${in_prefix}" \
      "${in_type}" \
      "${annotation}" \
      "${out_interval}" \
      "${out_prefix_csqs}" \
      "${chr}" )
  elif [ "${cluster}" = "sge" ]; then
    qsub -N "${ko_jname}" \
      -o "${ko_lname}.log" \
      -e "${ko_lname}.errors.log" \
      -P ${sge_project} \
      -wd $(pwd) \
      -t ${array_id} \
      -q "${sge_queue}" \
      -pe shmem ${nslots} \
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
  # note that this function waits for "ko_jname" (sge) which
  # is used as a dependency before starting job.
  submit_merge_job ${merge_prefix_regex} ${out_prefix_merge} ${ko_jname} ${slurm_jid}
}


submit_knockout_job "pLoF,damaging_missense" "1"




