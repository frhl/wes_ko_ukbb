#!/usr/bin/env bash
#
# author: Nbaya with revisions from flassen
#
#$ -N phase_chunks
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/phase_chunks.log
#$ -e logs/phase_chunks.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qe
#$ -t 20
#$ -V

set -o errexit
set -o nounset
module purge

# Set up
source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/phasing/04_phase_chunks.py"
readonly phasing_script="scripts/phasing/_phase_chunks.sh"
readonly spark_dir="data/tmp/spark"
set_up_hail
set_up_pythonpath_legacy

# Version of chrX-specific filter to use (options: females_only, both_sexes)
readonly chrX_filter_version="females_only"
# Number of variants within each interval
readonly min_interval_unit=1000
# Default size of phasing window in terms of variant count (should be a multiple of min_interval_unit)
readonly phasing_region_size=50000
# Minimum overlap between adjacent phasing windows
#readonly phasing_region_overlap=$(( ${phasing_region_size}/4 ))  
readonly phasing_region_overlap=10000 #$(( ${phasing_region_size}/4 ))  
# Maximum size of phasing window allowed, only used at the end of a chromosome
# Must be larger than phasing_region_size
readonly max_phasing_region_size=100000 

readonly chr=$( get_chr ${SGE_TASK_ID} )

# Cluster params
readonly software="shapeit4" #"shapeit4" or "eagle2"
readonly queue="short.qe"
readonly nslots=16

# what vcf should be phased
readonly vcf_dir=" data/unphased/wes_union_calls"
readonly vcf_to_phase="${vcf_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.bgz" # <--- change this back to eur!

# fam file for calculating switch errors
readonly pedigree_dir="/well/lindgren/UKBIOBANK/nbaya/resources"
readonly pedigree="${pedigree_dir}/ukb11867_pedigree.fam"

# Output paths
readonly out_dir="data/phased/wes_union_calls/chunks_non_eur"
readonly out_prefix="${out_dir}/ukb_non_eur_wes_union_calls_200k_chr${chr}"
readonly out_prefix_w_job_config="${out_prefix}-${nslots}x${queue}/${software}_prs${phasing_region_size}_pro${phasing_region_overlap}_mprs${max_phasing_region_size}"
readonly out="${out_prefix_w_job_config}.vcf.gz"
readonly out_symlink="${out_prefix}.vcf.gz"

readonly interval_dir="${out_dir}/intervals"
readonly interval_path="${interval_dir}/intervals_min_${min_interval_unit}_chr${chr}.tsv"
readonly phasing_interval_flags="--chrom ${chr} --min_interval_unit ${min_interval_unit}"

if [ -z "${interval_path}" ]; then
  raise_error "Getting intervals path failed"
fi

# Write phasing (minumum unit) intervals to slice later for phasing intervals
if [ ! -f ${interval_path} ]; then
  mkdir -p $( dirname ${interval_path} )
  SECONDS=0
  python3 ${hail_script} \
    ${phasing_interval_flags} \
    --write_intervals \
    --interval_path ${interval_path} \
    --target_vcf ${vcf_to_phase} \
    && print_update "Finished writing intervals for chr${chr}" ${SECONDS} \
    || raise_error "Writing intervals for chr${chr} failed" 
else
  print_update "${interval_path} already exists!"
fi

submit_phasing_job() {
  readonly max_phasing_idx=$( python3 ${hail_script} ${phasing_interval_flags} --phasing_region_size ${phasing_region_size} --phasing_region_overlap ${phasing_region_overlap} --max_phasing_region_size ${max_phasing_region_size} --get_max_phasing_idx --interval_path ${interval_path} )
  qsub -N "_c${chr}_${software}_phase_chunks" \
    -t 1-${max_phasing_idx} \
    -q ${queue} \
    -pe shmem ${nslots} \
    ${phasing_script} \
    ${chr} \
    ${vcf_to_phase} \
    ${min_interval_unit} \
    ${interval_path} \
    ${phasing_region_size} \
    ${phasing_region_overlap} \
    ${max_phasing_region_size} \
    ${out_prefix_w_job_config} \
    ${pedigree} \
    ${software}
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




