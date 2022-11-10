#!/usr/bin/env bash
#
# @description phase UK Biobank VCFs in optimal chunks
# @author Nbaya and flassen
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=phase_prephased_chunks
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/phase_prephased_chunks.log
#SBATCH --error=logs/phase_prephased_chunks.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=20

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly curwd=$(pwd)
readonly hail_script="scripts/phasing/09_phase_chunks.py"
readonly phasing_script="scripts/phasing/_phase_chunks.sh"
readonly spark_dir="data/tmp/spark"

# set +eu to avoid conda err
set +eu
set_up_hail
set_up_pythonpath_legacy
set -eu

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

# Version of chrX-specific filter to use (options: females_only, both_sexes)
readonly chrX_filter_version="females_only"
# Number of variants within each interval
readonly min_interval_unit=1000
# Default size of phasing window in terms of variant count (should be a multiple of min_interval_unit)
readonly phasing_region_size=100000
# Minimum overlap between adjacent phasing windows
readonly phasing_region_overlap=$(( ${phasing_region_size}/4 ))  
# Maximum size of phasing window allowed, only used at the end of a chromosome
# Must be larger than phasing_region_size
readonly max_phasing_region_size=150000
# minimum allele count allowed. Note that when using SHAPEIT4, singletons
# are randomly assigned an haplotype.
readonly pbwt_min_mac=2 # for shapeit5
readonly min_mac=2 # for shapeit4
# when prephased data is available, what is the phased set error?
# only available for shapeit4
readonly phased_set_error="0.0001" # 0.0001
# population effective size - only shapeit5.
readonly pop_effective_size=150000
# clsurm/sge parameters
readonly software="shapeit4" #"shapeit4", "shapeit5" or "eagle2"
readonly project="lindgren.prj"
readonly queue="short"
readonly nslots=16

# what vcf should be phased
readonly vcf_dir=" data/prephased/wes_union_calls"
readonly vcf_to_phase="${vcf_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.gz"

# SHAPEI5 requires a scaffold
readonly scaffold_dir="data/phased/calls/shapeit5/200k_from_500k"
readonly vcf_to_scaffold="${scaffold_dir}/NA"

# Output paths
readonly out_dir="data/phased/wes_union_calls/200k/whatshap/chunks/${software}"

readonly out_prefix="${out_dir}/ukb_wes_union_calls_whatshap_200k_chr${chr}"
readonly out_prefix_w_job_config="${out_prefix}-${nslots}x${queue}/${software}_prs${phasing_region_size}_pro${phasing_region_overlap}_mprs${max_phasing_region_size}"
readonly out="${out_prefix_w_job_config}.vcf.gz"
readonly out_symlink="${out_prefix}.vcf.gz"

readonly interval_dir="${out_dir}/intervals"
readonly interval_path="${interval_dir}/intervals_min_${min_interval_unit}_chr${chr}.tsv"
readonly phasing_interval_flags="--chrom ${chr} --min_interval_unit ${min_interval_unit}"

mkdir -p ${out_dir}

if [ -z "${interval_path}" ]; then
  raise_error "Getting intervals path failed"
fi

# Write phasing (minumum unit) intervals to slice later for phasing intervals
if [ ! -f ${interval_path} ]; then
  mkdir -p $( dirname ${interval_path} )
  SECONDS=0
  set -x
  python3 ${hail_script} \
    ${phasing_interval_flags} \
    --write_intervals \
    --interval_path ${interval_path} \
    --target_vcf ${vcf_to_phase} \
    && print_update "Finished writing intervals for chr${chr}" ${SECONDS} \
    || raise_error "Writing intervals for chr${chr} failed" 
  set +x
else
  print_update "${interval_path} already exists!"
fi

submit_phasing_job() {
  # get number of phasing indexes to run
  readonly max_phasing_idx=$( python3 ${hail_script} ${phasing_interval_flags} \
    --phasing_region_size ${phasing_region_size} \
    --phasing_region_overlap ${phasing_region_overlap} \
    --max_phasing_region_size ${max_phasing_region_size} \
    --get_max_phasing_idx \
    --interval_path ${interval_path} )
  # submit child script for phasing
  local slurm_tasks="1-${max_phasing_idx}"
  local slurm_jname="_c${chr}_${software}_phase_chunks"
  local slurm_lname="logs/_phase_chunks"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local slurm_nslots="${nslots}"
  set -x
  sbatch \
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
    ${phasing_script} \
    ${chr} \
    ${vcf_to_phase} \
    ${vcf_to_scaffold} \
    ${min_interval_unit} \
    ${interval_path} \
    ${phasing_region_size} \
    ${phasing_region_overlap} \
    ${max_phasing_region_size} \
    ${out_prefix_w_job_config} \
    ${software} \
    ${pbwt_min_mac} \
    ${min_mac} \
    ${phased_set_error} \
    ${pop_effective_size}
  set +x

}

if [ ! -f ${out} ]; then
  SECONDS=0
  module load BCFtools/1.12-GCC-10.3.0
  vcf_check ${vcf_to_phase}
  submit_phasing_job
  duration=${SECONDS}
  print_update "Finished submitting scattered phasing jobs for chr${chr}" "${duration}"
else
  print_update "Warning: ${out} already exists! Skipping." | tee /dev/stderr
fi




