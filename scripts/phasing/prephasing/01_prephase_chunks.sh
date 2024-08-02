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
#SBATCH --array=21

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly curwd=$(pwd)
readonly cluster=$( get_current_cluster)
readonly hail_script="scripts/phasing/prephasing/01_prephase_chunks.py"
readonly prephasing_script="scripts/phasing/prephasing/_prephase_chunks.sh"
readonly merge_script="scripts/phasing/prephasing/_prephase_merge.sh"
readonly spark_dir="data/tmp/spark"

# In each "chunk", we run whatshap across the samples and combine them.
readonly samples_per_chunk=100

# what matrix table should samples be drawn from. The only thing that
# this python-hail script is doing in creating a lists of samples in 
# chunks of 100s (as whatever parameter 'samples_per_chunk')
readonly input_samples_dir="data/prephased/wes_union_calls/intervals"
readonly input_samples_path="${input_samples_dir}/ukb_wes_union_calls_random_samples_50k_seed1995_chr21.mt/"
readonly input_samples_type="mt"

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

# Cluster params
# (note: errors after 50 min with 1000 samples/chunk with 4 slots)
readonly project="lindgren.prj"
readonly queue="short"
readonly nslots=2

# we need to provide whatshap the unphased VCF that we want to read-back phase.
readonly input_dir="data/unphased/wes_union_calls/prefilter/200k"
readonly input_path="${input_dir}/ukb_wes_union_calls_chr${chr}.vcf.gz" 
readonly input_type="vcf"

# Output paths
readonly out_dir="data/prephased/wes_union_calls/chunks/50k"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}"
readonly out_prefix_w_job_config="${out_prefix}_spc${samples_per_chunk}_${queue}/chr${chr}_chunk"
readonly out_merge_w_job_config="${out_prefix}_spc${samples_per_chunk}_${queue}.mergelist"

readonly out_merge_dir="data/prephased/wes_union_calls"
readonly out_merge_file="${out_merge_dir}/ukb_wes_union_calls_prephased_200k_chr${chr}.vcf.bgz"

# Interval paths
readonly interval_dir="${out_dir}/intervals"
readonly interval_path="${interval_dir}/intervals_spc${samples_per_chunk}_chr${chr}.tsv"

# Read locations (.cram)
readonly read_dir="/well/ukbb-wes/cram/oqfe/ukbb-11867"
readonly read_placeholder="${read_dir}/SAMPLE_oqfe.cram"

mkdir -p ${out_dir}
mkdir -p $( dirname ${out_prefix_w_job_config} )
mkdir -p $( dirname ${interval_path} )
mkdir -p ${out_merge_dir}

if [ ! -f ${interval_path} ]; then
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy
  python3 ${hail_script} \
    --write_interval \
    --samples_per_chunk ${samples_per_chunk} \
    --interval_path ${interval_path} \
    --input_path ${input_samples_path} \
    --input_type ${input_samples_type} \
    && print_update "Finished writing intervals for chr${chr}" ${SECONDS} \
    || raise_error "Writing intervals for chr${chr} failed" 
else
  echo "${interval_path} already exists!"
fi

get_max_interval_idx() {
  if [ -f ${interval_path} ]; then
    echo $( cat ${interval_path} | wc -l )
  else
    raise_error "${interval_path} does not exist."
  fi
}



submit_prephasing_job() {
  
  readonly max_interval_idx=$( get_max_interval_idx )
  readonly slurm_tasks="1-${max_interval_idx}"
  readonly slurm_jname="_c${chr}_prephase_chunks"
  readonly slurm_lname="logs/_prephase_chunks"
  readonly slurm_project="${project}"
  readonly slurm_queue="${queue}"
  readonly slurm_nslots="${nslots}"
  submit_prephasing_job_slurm

}

submit_prephasing_job_slurm() {
  echo "Submitting jobs with SLURM"
  readonly prephasing_jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname}.log" \
    --error="${slurm_lname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --array="${slurm_tasks}" \
    --open-mode="append" \
    --parsable \
    --constraint=skl-compat \
    ${prephasing_script} \
    ${chr} \
    ${input_path} \
    ${input_type} \
    ${interval_path} \
    ${max_interval_idx} \
    ${samples_per_chunk} \
    ${read_placeholder} \
    ${out_merge_w_job_config} \
    ${out_prefix_w_job_config} )
}

submit_merge_job() {
  local max_interval_idx=$( get_max_interval_idx )
  local slurm_jname="_c${chr}_prephase_merge"
  local slurm_lname="logs/_prephase_merge"
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
    --dependency="afterok:${prephasing_jid}" \
    --parsable \
    --constraint=skl-compat \
    ${merge_script} \
    ${out_prefix_w_job_config} \
    ${max_interval_idx} \
    ${out_merge_file}
  )
}


submit_prephasing_job "${cluster}"




