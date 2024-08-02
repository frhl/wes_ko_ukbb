#!/usr/bin/env bash
#
# @description phase UK Biobank VCFs in optimal chunks
# @author Nbaya and flassen
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=phase_chunks
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/phase_chunks.log
#SBATCH --error=logs/phase_chunks.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=21
#
#
#$ -N phase_chunks
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/phase_chunks.log
#$ -e logs/phase_chunks.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -t 21
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly curwd=$(pwd)
readonly hail_script="scripts/phasing/phasing/02_phase_chunks.py"
readonly phasing_script="scripts/phasing/phasing/_phase_chunks.sh"
readonly spark_dir="data/tmp/spark"

# set +eu to avoid conda err
set +eu
set_up_hail
set_up_pythonpath_legacy
set -eu

readonly cluster=$( get_current_cluster )
readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

## parmaters for chunks

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

# clsurm/sge parameters
readonly software="shapeit5" #"shapeit5" #"shapeit4" or "eagle2"
readonly project="lindgren.prj"
readonly queue="long"
readonly nslots=16

## paramters for phasing with shapeit
readonly phased_set_error="0.0001" # 0.0001
readonly pbwt_min_mac=2 # for shapeit5r
readonly pbwt_depth=2 # 5
readonly pbwt_modulo=0.1 # default is 0.1 but 0.0004 ( 0.2 / 50 ) is default value when using --sequencing arugment
readonly pbwt_mdr=0.1
readonly pop_effective_size=15000

readonly tranche="200k"

# what vcf should be phased
readonly vcf_dir="data/unphased/wes_union_calls/prefilter/${tranche}"
readonly vcf_to_phase="${vcf_dir}/ukb_wes_union_calls_chr${chr}.vcf.gz"

# SHAPEI5 requires a scaffold
readonly scaffold_dir="data/phased/wes_union_calls/${tranche}/shapeit5/phase_common/newrun"
readonly vcf_to_scaffold="${scaffold_dir}/ukb_wes_union_calls_${tranche}_chr${chr}_phase_common.vcf.gz"

#readonly scaffold_dir="data/phased/calls/shapeit5/200k_from_500k"
#readonly vcf_to_scaffold="${scaffold_dir}/ukb_phased_calls_200k_from_500k_chr${chr}.vcf.bgz"

# Output paths
#readonly out_dir="data/phased/wes_union_calls/200k/eagle2/shapeit4/chunks"
readonly out_dir="data/phased/wes_union_calls/${tranche}/shapeit5/phase_rare_test"

readonly out_prefix="${out_dir}/ukb_wes_union_calls_${software}_${tranche}_chr${chr}"
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
  if [ "${cluster}" == "slurm" ]; then
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
      ${pbwt_depth} \
      ${pbwt_modulo} \
      ${pbwt_mdr} \
      ${phased_set_error} \
      ${pop_effective_size}
    set +x
  elif [ "${cluster}" == "sge" ]; then
    qsub -N "${slurm_jname}" \
      -o "${slurm_lname}.log" \
      -e "${slurm_lname}.errors.log" \
      -t ${slurm_tasks} \
      -q "short.qc" \
      -pe shmem ${slurm_nslots} \
      -wd $(pwd) \
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
      ${pbwt_depth} \
      ${pbwt_modulo} \
      ${pbwt_mdr} \
      ${phased_set_error} \
      ${pop_effective_size}
  else
    >&2 echo "${cluster} is not valid!"
  fi

}

if [ ! -f ${out} ]; then
  SECONDS=0
  module load BCFtools/1.12-GCC-10.3.0
  vcf_check ${vcf_to_phase}
  vcf_check ${vcf_to_scaffold}
  make_tabix ${vcf_to_phase}
  make_tabix ${vcf_to_scaffold}
  submit_phasing_job
  duration=${SECONDS}
  print_update "Finished submitting scattered phasing jobs for chr${chr}" "${duration}"
else
  print_update "Warning: ${out} already exists! Skipping." | tee /dev/stderr
fi




