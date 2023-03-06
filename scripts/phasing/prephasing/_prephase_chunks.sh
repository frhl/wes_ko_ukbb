#!/usr/bin/env bash
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly cluster=$( get_current_cluster)
readonly project="lindgren.prj"
readonly queue="short"
readonly nslots="1"

readonly chr=${1?Error: Missing arg1 (chr)} # Chromosome, e.g. "1" for chrom 1
readonly input_path=${2?Error: Missing arg2 (input_path)} 
readonly input_type=${3?Error: Missing arg3 (input_type)} 
readonly interval_path=${4?Error: Missing arg4 (intervals_path)} 
readonly max_interval_idx=${5?Error: Missing arg5 (intervals_path)} 
readonly samples_per_chunk=${6?Error: Missing arg5 (intervals_path)} 
readonly read_placeholder=${7?Error: Missing arg5 (intervals_path)} 
readonly main_merge_file=${8?Error: Missing arg6 ()} 
readonly out_prefix=${9?Error: Missing arg6 ()} 

readonly prephase_sample_script="scripts/phasing/prephasing/_prephase_sample.sh"
readonly hail_script="scripts/phasing/prephasing/01_prephase_chunks.py"
readonly merge_script="scripts/phasing/prephasing/_prephase_merge.sh"

readonly chunk_idx=$( get_array_task_id ) # one-based index for which phasing interval to phase
readonly out_prefix_w_interval_idx="${out_prefix}_${chunk_idx}of${max_interval_idx}"
readonly splitted="${out_prefix_w_interval_idx}"
readonly splitted_input="${splitted}.vcf.bgz"
readonly splitted_type="vcf"
# parameters for merge
readonly merge_list="${out_prefix_w_interval_idx}.mergelist"
readonly merge_type="vcf"
readonly out_merge_file="${out_prefix_w_interval_idx}_prephased"
readonly out_merge_type="vcf"

# append files to be merged
echo "main mergefile: ${main_merge_file}"
echo "${out_merge_file}.vcf.gz" >> ${main_merge_file}


#rm -f ${merge_list}
mkdir -p ${out_prefix_w_interval_idx}

split_to_chunks() {
  # use hail to split to pre-defined chunks of samples
  SECONDS=0
  local current_interval=${1}
  local file_to_split=${2}
  local type_to_split=${3}
  local out_file=${4}
  local out_type="vcf" # always vcf
  if [ ! -f "${out_file}.vcf.bgz" ]; then
    module purge
    set_up_hail
    set_up_pythonpath_legacy
    python3 ${hail_script} \
      --input_path ${file_to_split} \
      --input_type ${type_to_split} \
      --output_path ${out_file} \
      --output_type ${out_type} \
      --interval_path ${interval_path} \
      --interval_idx ${current_interval} \
      --split_by_interval \
      && print_update "Finished splitting for ${out_prefix_w_interval_idx}" ${SECONDS} \
      || raise_error "Error when splitting chunks for ${out_prefix_w_interval_idx}"
    set +x
  else
    >&2 echo "${out_file} (split) already exists. Skipping!"
  fi
}

submit_prephasing_sample_job() {

  local slurm_tasks="1-${samples_per_chunk}"
  local slurm_jname="_c${chr}_i${chunk_idx}"
  local slurm_lname="logs/_prephase_sample"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local slurm_nslots="${nslots}"
  if [ "${cluster}" == "slurm" ]; then
    submit_prephasing_job_slurm
  elif [ "${cluster}" == "sge" ]; then
    submit_prephasing_job_sge
  else
    echo "${cluster} is not valid!"
  fi

}


submit_prephasing_job_slurm() {
  echo "Submitting jobs with SLURM"
  readonly prephasing_jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname}.log" \
    --error="${slurm_lname}.errors.log" \
    --chdir="$(pwd)" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --array="${slurm_tasks}" \
    --open-mode="append" \
    --parsable \
    --constraint=skl-compat \
    ${prephase_sample_script} \
    ${splitted_input} \
    ${interval_path} \
    ${chunk_idx} \
    ${read_placeholder} \
    ${merge_list} \
    ${out_prefix_w_interval_idx} )
}

submit_prephasing_job_sge() {
  echo "Submitting jobs with SGE"
  qsub -N "${slurm_jname}" \
    -o "${slurm_lname}.log" \
    -e "${slurm_lname}.errors.log" \
    -t ${slurm_tasks} \
    -q "short.qc" \
    -pe shmem ${slurm_nslots} \
    -wd $(pwd) \
    ${prephase_sample_script} \
    ${splitted_input} \
    ${interval_path} \
    ${chunk_idx} \
    ${read_placeholder} \
    ${merge_list} \
    ${out_prefix_w_interval_idx}
}

submit_merge_job_slurm() {
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
    --chdir="$(pwd)" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --parsable \
    --constraint=skl-compat \
    ${merge_script} \
    ${merge_list} \
    ${merge_type} \
    ${out_merge_file} \
    ${out_merge_type} 
  )
  #--dependency="afterok:${prephasing_jid}" \
}

submit_merge_job_sge() {
  local slurm_jname="_c${chr}_prephase_merge"
  local slurm_lname="logs/_prephase_merge"
  local wait_for="_c${chr}_i${chunk_idx}"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local slurm_nslots="1"
  qsub -N "${slurm_jname}" \
    -o "${slurm_lname}.log" \
    -e "${slurm_lname}.errors.log" \
    -q "short.qc" \
    -pe shmem ${slurm_nslots} \
    -hold_jid "${wait_for}" \
    -wd $(pwd) \
    ${merge_script} \
    ${merge_list} \
    ${merge_type} \
    ${out_merge_file} \
    ${out_merge_type} 
}

# check if phasing and merging has already been completed
force_rm_bad_vcf "${out_merge_file}.vcf.gz" 
if [ ! -f "${out_merge_file}.vcf.gz" ]; then
  
  # split main MatrixTable into 
  # chunks of equally sized VCFs by sample
  if [ ! -f "${splitted_input}" ]; then
    split_to_chunks \
      ${chunk_idx} \
      ${input_path} \
      ${input_type} \
      ${splitted}
  fi

  # make tabix of VCF
  if [ ! -f "${splitted_input}.tbi" ]; then
    module purge
    module load BCFtools/1.12-GCC-10.3.0
    make_tabix "${splitted_input}" "tbi"
  fi

  # only submit jobs if the right files exists
  if [ -f "${splitted_input}.tbi" ]; then
    # submit workers for each sample
    #submit_prephasing_sample_job ${cluster}
    submit_merge_job_slurm
  else
    >&2 echo "Error: Prephasing job was not submitted because VCF/VCF.tbi files do not exists."
  fi


else
  >&2 echo "${out_merge_file}.vcf.gz already exists. Skipping."
fi

