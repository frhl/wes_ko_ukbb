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

readonly prephase_sample_script="scripts/phasing/_prephase_sample.sh"
readonly hail_script="scripts/phasing/05_prephase_chunks.py"
readonly rscript="scripts/phasing/_prephase_sample.R"

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


extract_sample() {
    local in_vcf=${1}
    local sample=${2}
    local out_vcf=${3}
    if [ ! -f "${out_vcf}.gz" ]; then
      module load BCFtools/1.12-GCC-10.3.0
      bcftools view -s ${sample} -o ${out_vcf} ${in_vcf} && bgzip ${out_vcf}
      make_tabix "${out_vcf}.gz" "tbi"
    fi
}

prephase_sample() {

  local eid=${1}
  local to_phase_vcf=${2}
  local read_path=${read_placeholder/SAMPLE/${eid}}
  local path_unphased="${out_prefix}/s${sample_idx}_eid${eid}.vcf"
  local path_unphased_gz="${path_unphased}.gz"
  local prefix_phased="${out_prefix}/s${sample_idx}_eid${eid}_phased"
  local path_phased_vcf="${prefix_phased}.vcf"
  local path_phased_gz="${prefix_phased}.vcf.gz"
  local out_reads="${out_prefix}/s${sample_idx}_eid${eid}.reads"
  local tmp_reads="${out_prefix}/s${sample_idx}_eid${eid}.reads.tmp"

  # get single sample
  extract_sample \
    ${to_phase_vcf} \
    ${eid} \
    ${path_unphased}

  # prephase single sample
  SECONDS=0
  whatshap phase \
    --reference="${grch38}" \
    --output="${path_phased_vcf}" \
    --output-read-list="${out_reads}" \
    --ignore-read-groups \
    --indels \
    ${path_unphased_gz} \
    ${read_path} \
    && print_update "Finished prephasing ${path_unphased}" ${SECONDS} \
    || raise_error "Error prephasing ${path_unphased}"
  
  # bgzip sample
  module load BCFtools/1.12-GCC-10.3.0
  bgzip ${path_phased_vcf}
  make_tabix "${path_phased_gz}" "tbi"

}



# check if phasing and merging has already been completed
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
 
  # need modules for VCF handling and
  # the module which has whatshap installed
  module purge
  module load HTSlib/1.12-GCC-10.3.0
  module load BCFtools/1.12-GCC-10.3.0
  set_up_whatshap
 
  for sample_idx in $( seq 1 ${samples_per_chunk}); do 
    eid=$( Rscript ${rscript} --chunk_idx ${chunk_idx} --sample_idx ${sample_idx} --interval_path ${interval_path} )
    prephase_sample ${eid} ${splitted_input} 
  done


else
  >&2 echo "${out_merge_file}.vcf.gz already exists. Skipping."
fi

