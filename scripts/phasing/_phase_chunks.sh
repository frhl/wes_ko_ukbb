#!/usr/bin/env bash
#
#$ -N _phase_chunks
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_phase_chunks.log
#$ -e logs/_phase_chunks.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qa
#$ -V

set -o errexit
set -o nounset
module purge

# Set up

#source /well/lindgren/UKBIOBANK/nbaya/resources/ukb_utils/bash/qsub_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh
source_ukb_utils_scripts hail vcf phasing
spark_dir="data/tmp/spark"
set_up_hail
set_up_pythonpath_legacy
#add_module_to_pythonpath ukb_utils ukb_wes_qc phase_ukb_imputed

# args
readonly chr=${1?Error: Missing arg1 (chr)} # Chromosome, e.g. "1" for chrom 1
readonly vcf_to_phase=${2?Error: Missing arg2 (vcf_to_phase)} # Path of VCF to phase
readonly min_interval_unit=${3?Error: Missing arg3 (min_interval_unit)} 
readonly phasing_region_size=${4?Error: Missing arg4 (phasing_region_size)} # Minimum guaranteed size of phasing window in terms of variant count
readonly phasing_region_overlap=${5?Error: Missing arg5 (phasing_region_overlap)} # Minimum overlap between adjacent phasing windows
readonly max_phasing_region_size=${6?Error: Missing arg6 (max_phasing_region_size)} # Maximum size of phasing window allowed, only used at the end of a chromosome. Must be larger than phasing_region_size
readonly out_prefix=${7?Error: Missing arg7 (path prefix for output intermediate VCF)} # Path to output phased VCF
#readonly scaffold=${8?Error: Missing arg8 (scaffold)} # Path of VCF to use as haplotype scaffold

#readonly shapeit_version="4.2.2"
readonly interval_flags="--chrom ${chr} --min_interval_unit ${min_interval_unit} --phasing_region_size ${phasing_region_size} --phasing_region_overlap ${phasing_region_overlap} --max_phasing_region_size ${max_phasing_region_size}"

readonly hail_script="scripts/phasing/_phase_chunks.py"

readonly phasing_idx=${SGE_TASK_ID} # one-based index for which phasing interval to phase

readonly max_phasing_idx=$( python3 ${hail_script} ${interval_flags} --get_max_phasing_idx )
readonly out_prefix_w_phasing_idx="${out_prefix}.${phasing_idx}of${max_phasing_idx}"
readonly out="${out_prefix_w_phasing_idx}.vcf.gz"
readonly log="${out_prefix_w_phasing_idx}.log"


phase_with_shapeit() {
  SECONDS=0
  local gmap="/well/lindgren/flassen/software/SHAPEIT4/b38.gmap/chr${chr}.b38.gmap.gz"
  mkdir -p $( dirname ${out} )
  readonly region=$( python3 ${hail_script} ${interval_flags} --get_interval --phasing_idx ${phasing_idx} )
  print_update "Starting phasing for ${region}, out: ${out}"
  set -x
  shapeit4.2 \
    --input ${vcf_to_phase} \
    --map ${gmap} \
    --region "chr${region}" \
    --thread $(( ${NSLOTS}-1 )) \
    --sequencing \
    --log ${log} \
    --output ${out}
  set +x
  duration=${SECONDS}
  print_update "Finished phasing non-singleton variants for ${region} (chr${chr} ${phasing_idx}/${max_phasing_idx}, out: ${out}" "${duration}"
}

if [ ! -f ${out} ]; then
  
  module load SHAPEIT4/4.2.2-foss-2021a
  phase_with_shapeit  
else
  print_update "Warning: ${out} already exists! Skipping." | tee /dev/stderr
fi

make_tabix ${out}



