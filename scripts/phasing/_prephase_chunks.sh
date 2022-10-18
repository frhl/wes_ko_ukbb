#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"

readonly chr=${1?Error: Missing arg1 (chr)} # Chromosome, e.g. "1" for chrom 1
readonly input_path=${2?Error: Missing arg2 (input_path)} 
readonly input_type=${3?Error: Missing arg3 (input_type)} 
readonly interval_path=${4?Error: Missing arg4 (intervals_path)} 
readonly max_interval_idx=${5?Error: Missing arg5 (intervals_path)} 
readonly read_placeholder=${6?Error: Missing arg5 (intervals_path)} 
readonly out_prefix=${7?Error: Missing arg6 ()} 

readonly hail_script="scripts/phasing/03_prephase_chunks.py"
readonly rscript="scripts/phasing/_prephase_chunks.R"

readonly interval_idx=${SLURM_ARRAY_TASK_ID} # one-based index for which phasing interval to phase
readonly out_prefix_w_interval_idx="${out_prefix}.${interval_idx}of${max_interval_idx}"

readonly out_splitted="${out_prefix_w_interval_idx}"
readonly out_splitted_vcf="${out_prefix_w_interval_idx}.vcf.bgz"
readonly out_phased="${out_prefix_w_interval_idx}.phased.vcf.gz"
readonly log="${out_prefix_w_interval_idx}.log"

# setup reference (Required for whatshap) 
readonly refdir="/well/lindgren/flassen/ressources/genome_reference/broad"
readonly grch38="${refdir}/Homo_sapiens_assembly38.fasta"
export REF_PATH="${refdir}/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
export REF_CACHE="${refdir}/ref/cache/%2s/%2s/%s"

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
  else
    >&2 echo "${out_file} (split) already exists. Skipping!"
  fi
}

get_read_files() {
  # use R to get paths to relevant read files
  module purge
  set_up_rpy
  local current_interval=${1}
  local read_sample_path=${2}
  local reads=$( Rscript ${rscript} \
    --interval_idx ${current_interval} \
    --interval_path ${interval_path} \
    --read_placeholder ${read_sample_path} )
  echo ${reads}
}

prephase_with_whatshap() {
  # perform prephasing using whatshap
  SECONDS=0
  local vcf_to_phase=${1}
  local vcf_result=${2}
  local cram_files=${3}
  module purge
  module load htslib/1.8-gcc5.4.0
  set_up_whatshap
  set -x
  whatshap phase \
    --reference="${grch38}" \
    --output="${vcf_result}" \
    --indels \
    ${vcf_to_phase} \
    ${cram_files} \
    && print_update "Finished prephasing ${out_prefix_w_interval_idx}" ${SECONDS} \
    || raise_error "Error prephasing ${out_prefix_w_interval_idx}"
  set +x
}



if [ ! -f ${out_phased} ]; then
  
  # split main MatrixTable into 
  # chunks of equally sized VCFs by sample
  split_to_chunks \
    ${interval_idx} \
    ${input_path} \
    ${input_type} \
    ${out_splitted}

  # make tabix of VCF
  if [ ! -f "${out_splitted_vcf}.tbi" ]; then
    module purge
    module load BCFtools/1.12-GCC-10.3.0
    make_tabix "${out_splitted_vcf}" "tbi"
  fi
  
  # get variable with paths to relevant read files
  readonly reads=$( get_read_files ${interval_idx} ${read_placeholder} )
  
  # perform read-backed phasing using whatshap
  prephase_with_whatshap \
    ${out_splitted_vcf} \
    ${out_phased} \
    ${reads}
    

else
  print_update "Warning: ${out_phased} already exists! Skipping." | tee /dev/stderr
fi




