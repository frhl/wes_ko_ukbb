#!/usr/bin/env bash
#
# Nbaya
#
#$ -N phase_chunks
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/phase_chunks.log
#$ -e logs/phase_chunks.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qe
#$ -t 21
#$ -V

set -o errexit
set -o nounset
module purge

# Set up
source utils/qsub_utils.sh
#source /well/lindgren/UKBIOBANK/nbaya/resources/ukb_utils/bash/qsub_utils.sh
source utils/hail_utils.sh
#source_ukb_utils_scripts hail vcf

spark_dir="data/tmp/spark"
set_up_hail
set_up_pythonpath_legacy
#add_module_to_pythonpath ukb_utils ukb_wes_qc phase_ukb_imputed phase_ukb_wes

readonly hail_script="scripts/phasing/_phase_chunks.py"
readonly utils_script="scripts/phasing/_phase_chunks_utils.py"

readonly chr=$( get_chr ${SGE_TASK_ID} )

# Version of chrX-specific filter to use (options: females_only, both_sexes)
readonly chrX_filter_version="females_only"
# Number of variants within each interval
readonly min_interval_unit=1000
# Default size of phasing window in terms of variant count (should be a multiple of min_interval_unit)
readonly phasing_region_size=50000 #100000
# Minimum overlap between adjacent phasing windows
readonly phasing_region_overlap=$(( ${phasing_region_size}/4 ))  
# Maximum size of phasing window allowed, only used at the end of a chromosome
# Must be larger than phasing_region_size
readonly max_phasing_region_size=100000 

readonly phasing_interval_flags="--chrom ${chr} --min_interval_unit ${min_interval_unit}"
readonly intervals_path=$( python3 ${utils_script} ${phasing_interval_flags} --print_phasing_intervals_path )

if [ -z "${intervals_path}" ]; then
  raise_error "Getting intervals path failed"
fi


# Write phasing (minumum unit) intervals to slice later for phasing intervals
if [ ! -f ${intervals_path} ]; then
  mkdir -p $( dirname ${intervals_path} )
  set_up_hail
  set_up_pythonpath_legacy
  SECONDS=0
  python3 ${hail_script} \
    ${phasing_interval_flags} \
    --write_intervals \
    && print_update "Finished writing intervals for chr${chr}" ${SECONDS} \
    || raise_error "Writing intervals for chr${chr} failed" 
else
  print_update "${intervals_path} already exists!"
fi


# Cluster params
readonly queue="short.qa"
readonly nslots=10

# Phased VCF to use as a scaffold for SHAPEIT4
readonly vcf_dir=" data/unphased/wes_union_calls"
readonly vcf_to_phase="${vcf_dir}/ukb_eur_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly scaffold=""

# Output paths
readonly out_dir="data/phased/wes/call_union"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}"
readonly out_prefix_w_job_config="${out_prefix}-${nslots}x${queue}/prs${phasing_region_size}_pro${phasing_region_overlap}_mprs${max_phasing_region_size}"
readonly out="${out_prefix_w_job_config}.vcf.gz"
readonly out_symlink="${out_prefix}.vcf.gz"

# Phasing script
readonly phasing_script="scripts/phasing/_phase_chunks.sh"

# Add BCFtools to PATH
export PATH="${PATH}:/gpfs3/users/gms/whv244/bcftools_1.11/bin"
#bcftools_check


submit_phasing_job() {
  readonly max_phasing_idx=$( python3 ${hail_script} ${phasing_interval_flags} --phasing_region_size ${phasing_region_size} --phasing_region_overlap ${phasing_region_overlap} --max_phasing_region_size ${max_phasing_region_size} --get_max_phasing_idx )
  qsub -N "_c${chr}_phase_chunks" \
    -t 1-${max_phasing_idx} \
    -q ${queue} \
    -pe shmem ${nslots} \
    ${phasing_script} \
    ${chr} \
    ${vcf_to_phase} \
    ${min_interval_unit} \
    ${phasing_region_size} \
    ${phasing_region_overlap} \
    ${max_phasing_region_size} \
    ${out_prefix_w_job_config} \
    ${scaffold}
}

if [ ! -f ${out} ]; then
  SECONDS=0

  #vcf_check ${vcf_to_phase}
  #vcf_check ${scaffold}
  submit_phasing_job

  duration=${SECONDS}
  print_update "Finished submitting scattered phasing jobs for chr${chr}" "${duration}"
else
  print_update "Warning: ${out} already exists! Skipping." | tee /dev/stderr
fi




