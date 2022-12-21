#!/usr/bin/env bash
#
# @description: amalgamate variants by phase to infer knockouts by genes.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=encode_vcf_array
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/encode_vcf_array.log
#SBATCH --error=logs/encode_vcf_array.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-2
#
#
#$ -N encode_vcf_array
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/encode_vcf_array.log
#$ -e logs/encode_vcf_array.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1-22
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly curwd=$(pwd)
readonly spark_dir="data/tmp/spark"
readonly bash_script="scripts/_encode_vcf_array.sh"
readonly merge_script="scripts/_merge_encode_vcf_array.sh"
readonly hail_script="scripts/_write_gene_intervals.py"

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_dir="data/mt/prefilter/pp90"
readonly merge_dir="data/knockouts/alt/pp90/extracted_array"
readonly out_dir="data/knockouts/alt/pp90/extracted_array/chr${chr}"
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
readonly slurm_queue="short"
readonly sge_queue="short.qc" 

mkdir -p ${out_dir}
mkdir -p ${merge_dir}

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

# How many chunks should be run
extract_genes
readonly genes_per_chunk=100 
readonly n_genes=$( cat ${out_interval} | wc -l )
readonly chunks=$(( (${genes_per_chunk}+${n_genes}-1) / ${genes_per_chunk} ))
readonly array_id="1-${chunks}"


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
      "${chunks}" \
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
      "${chunks}" \
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
  local merge_prefix_regex="${out_prefix_csqs}_ESNG"

  # slurm specific paramters 
  local ko_jname="_c${chr}_extr_${annotation}"
  local ko_lname="logs/_encode_vcf_array"
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
      "${genes_per_chunk}" \
      "${chunks}" \
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
      "${genes_per_chunk}" \
      "${chunks}" \
      "${chr}"
  else
    >&2 echo "${cluster} is not a valid cluster."
  fi
  # note that this function waits for "ko_jname" (sge) which
  # is used as a dependency before starting job.
  submit_merge_job ${merge_prefix_regex} ${out_prefix_merge} ${ko_jname} ${slurm_jid}
}


#submit_knockout_job "pLoF,damaging_missense" "2"
#submit_knockout_job "damaging_missense" "2"
#submit_knockout_job "pLoF" "2"
submit_knockout_job "synonymous" "2"
#submit_knockout_job "other_missense" "3"




