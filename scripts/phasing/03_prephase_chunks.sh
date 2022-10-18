#!/usr/bin/env bash
#
# @description split into chunks of samples that are then pre-phased using whatshap
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=prephase_chunks
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/prephase_chunks.log
#SBATCH --error=logs/prephase_chunks.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=22

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly curwd=$(pwd)
readonly hail_script="scripts/phasing/03_prephase_chunks.py"
readonly prephasing_script="scripts/phasing/_prephase_chunks.sh"
readonly merge_script="scripts/phasing/_prephase_merge.sh"
readonly spark_dir="data/tmp/spark"

# how many samples should there be in each chunk 
readonly samples_per_chunk=1000
readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} )

# Cluster params
readonly project="lindgren.prj"
readonly queue="short"
readonly nslots=4

# what file should be split up
readonly input_dir=" data/unphased/wes_union_calls"
readonly input_path="${input_dir}/ukb_wes_union_calls_200k_chr${chr}.mt" 
readonly input_type="mt"

# Output paths
readonly out_dir="data/phased/wes_union_calls/prephased/chunks"
readonly out_prefix="${out_dir}/ukb_eur_wes_union_calls_200k_chr${chr}"
readonly out_prefix_w_job_config="${out_prefix}-${queue}/chr${chr}_spc${samples_per_chunk}"

readonly out_merge_dir="data/phased/wes_union_calls/prephased"
readonly out_merge_file="${out_merge_dir}/ukb_eur_wes_union_calls_prephased_200k_chr${chr}.vcf.bgz"

# Interval paths
readonly interval_dir="${out_dir}/intervals"
readonly interval_path="${interval_dir}/intervals_spc${samples_per_chunk}_chr${chr}.tsv"

# Read locations
readonly read_dir="/well/ukbb-wes/cram/oqfe/ukbb-11867"
readonly read_placeholder="${read_dir}/SAMPLE_oqfe.cram"

mkdir -p ${out_dir}


if [ ! -f ${interval_path} ]; then
  mkdir -p $( dirname ${interval_path} )
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy
  python3 ${hail_script} \
    --write_interval \
    --samples_per_chunk ${samples_per_chunk} \
    --interval_path ${interval_path} \
    --input_path ${input_path} \
    --input_type ${input_type} \
    && print_update "Finished writing intervals for chr${chr}" ${SECONDS} \
    || raise_error "Writing intervals for chr${chr} failed" 
else
  print_update "${interval_path} already exists!"
fi

get_max_interval_idx() {
  if [ -f ${interval_path} ]; then
    echo $( cat ${interval_path} | wc -l )
  else
    raise_error "${interval_path} does not exist."
  fi
}



submit_prephasing_job() {
  # get number of phasing indexes to run
  local max_interval_idx=$( get_max_interval_idx )
  local slurm_tasks="1-2" #-${max_phasing_idx}"
  local slurm_jname="_c${chr}_prephase_chunks"
  local slurm_lname="logs/_prephase_chunks"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local slurm_nslots="${nslots}"
  local jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname}.log" \
    --error="${slurm_lname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --array="${slurm_tasks}" \
    --open-mode="append" \
    --constraint=skl-compat \
    ${prephasing_script} \
    ${chr} \
    ${input_path} \
    ${input_type} \
    ${interval_path} \
    ${max_interval_idx} \
    ${read_placeholder} \
    ${out_prefix_w_job_config} )
  echo ${jid}
}

submit_merge_job() {
  # merge resulting chunk files
  local dependency=${1}
  local max_interval_idx=$( get_max_interval_idx )
  local slurm_jname="_c${chr}_merge_prephased_chunks"
  local slurm_lname="logs/_merge_prephased_chunks"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local slurm_nslots="2"
  local jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname}.log" \
    --error="${slurm_lname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --dependency="afterok:${dependency}" \
    --constraint=skl-compat \
    ${merge_script} \
    ${out_prefix_w_job_config} \
    ${max_interval_idx} \
    ${out_merged_file}
  )

}


readonly prephasing_jid=$( submit_prephasing_job )




