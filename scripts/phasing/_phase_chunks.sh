#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
set_up_hail
set_up_pythonpath_legacy

readonly chr=${1?Error: Missing arg1 (chr)} # Chromosome, e.g. "1" for chrom 1
readonly vcf_to_phase=${2?Error: Missing arg2 (vcf_to_phase)} # Path of VCF to phase
readonly vcf_to_scaffold=${3?Error: Missing arg2 (vcf_to_phase)} # Path of VCF to phase
readonly min_interval_unit=${4?Error: Missing arg3 (min_interval_unit)} 
readonly interval_path=${5?Error: Missing arg4 (min_interval_unit)} 
readonly phasing_region_size=${6?Error: Missing arg5 (phasing_region_size)} # Minimum guaranteed size of phasing window in terms of variant count
readonly phasing_region_overlap=${7?Error: Missing arg6 (phasing_region_overlap)} # Minimum overlap between adjacent phasing windows
readonly max_phasing_region_size=${8?Error: Missing arg7 (max_phasing_region_size)} # Maximum size of phasing window allowed, only used at the end of a chromosome. Must be larger than phasing_region_size
readonly out_prefix=${9?Error: Missing arg8 (path prefix for output intermediate VCF)} # Path to output phased VCF
readonly software=${10?Error: Missing arg10 (software)} # Path of VCF to use as haplotype scaffold
readonly threads=$(( ${SLURM_CPUS_ON_NODE} - 1))

readonly hail_script="scripts/phasing/09_phase_chunks.py"
readonly interval_flags="--chrom ${chr} --min_interval_unit ${min_interval_unit} --phasing_region_size ${phasing_region_size} --phasing_region_overlap ${phasing_region_overlap} --max_phasing_region_size ${max_phasing_region_size}"
readonly phasing_idx=$( get_array_task_id ) # one-based index for which phasing interval to phase

readonly max_phasing_idx=$( python3 ${hail_script} ${interval_flags} --get_max_phasing_idx --interval_path ${interval_path})
readonly out_prefix_w_phasing_idx="${out_prefix}.${phasing_idx}of${max_phasing_idx}"
readonly out="${out_prefix_w_phasing_idx}.vcf.gz"
readonly log="${out_prefix_w_phasing_idx}.log"

mkdir -p $( dirname ${out} )

phase_with_shapeit4() {
  SECONDS=0
  readonly min_mac=2 # default value. Unable to phase singletons.
  readonly ps_error_rate=0.0001 # default value
  readonly gmap="/well/lindgren/flassen/software/SHAPEIT4/b38.gmap/chr${chr}.b38.gmap.gz"
  readonly region=$( python3 ${hail_script} ${interval_flags} --get_interval --phasing_idx ${phasing_idx} --interval_path ${interval_path} )
  print_update "Starting SHAPEIT5 phasing for ${region} (min_mac >= ${min_mac}) out: ${out}"
  set -x
  set_up_shapeit5
  ${SHAPEIT_phase_rare} \
    --input-plain ${vcf_to_phase} \
    --scaffold ${vcf_to_scaffold} \
    --input-region "chr${region}" \
    --scaffold-region "chr${region}" \
    --pbwt-mac ${min_mac} \
    --map ${gmap} \
    --thread ${threads} \
    --log ${log} \
    --output ${out}
  set +x
  duration=${SECONDS}
  print_update "Finished phasing variants for ${region} (chr${chr} ${phasing_idx}/${max_phasing_idx}, out: ${out}" "${duration}"
}

phase_with_shapeit5() {
  SECONDS=0
  readonly min_mac=2 # default value. Unable to phase singletons.
  readonly gmap="/well/lindgren/flassen/software/SHAPEIT4/b38.gmap/chr${chr}.b38.gmap.gz"
  readonly region=$( python3 ${hail_script} ${interval_flags} --get_interval --phasing_idx ${phasing_idx} --interval_path ${interval_path} )
  print_update "Starting SHAPEIT5 phasing for ${region} (min_mac >= ${min_mac}) out: ${out}"
  set -x
  ${SHAPEIT_phase_common} \
    --input ${vcf_to_phase} \
    --map ${gmap} \
    --region "chr${region}" \
    --thread ${threads} \
    --pbwt-mac ${min_mac} \
    --output ${out} \
    --log ${log} \
    && print_update "Finished phasing variants for chr${chr}, out: ${out}" "${SECONDS}" \
    || raise_error "$( print_update "Phasing variants failed for chr${chr}" ${SECONDS} )" 
  set +x
  duration=${SECONDS}
  print_update "Finished phasing variants for ${region} (chr${chr} ${phasing_idx}/${max_phasing_idx}, out: ${out}" "${duration}"
}


phase_with_eagle2() {
  SECONDS=0
  readonly eagle="/well/lindgren/flassen/software/eagle/Eagle_v2.4.1/./eagle"
  readonly gmap="/well/lindgren/flassen/software/eagle/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
  readonly region=$( python3 ${hail_script} ${interval_flags} --get_interval --phasing_idx ${phasing_idx} --interval_path ${interval_path} )
  readonly bp_start="$( echo "${region%-*}" | sed 's/^.*://' )"
  readonly bp_end="${region##*-}"
  print_update "Starting eagle2 phasing in region ${bp_start}-${bp_end}, out: ${out}"
  set -x
  ${eagle} \
    --vcf ${vcf_to_phase} \
    --bpStart ${bp_start} \
    --bpEnd ${bp_end} \
    --outPrefix ${out_prefix_w_phasing_idx} \
    --geneticMapFile ${gmap} \
    --numThreads ${threads} \
    --maxMissingPerIndiv=0.1 \
    --maxMissingPerSnp=0.1 \
    --Kpbwt=20000 \
    --pbwtOnly
  set +x
  duration=${SECONDS}
  print_update "Finished phasing for ${region} (chr${chr} ${phasing_idx}/${max_phasing_idx}, out: ${out}" "${duration}"
}



if [ ! -f ${out} ]; then
  if [ ${software} = "shapeit4" ]; then
    module load SHAPEIT4/4.2.2-foss-2021a
    phase_with_shapeit4
  elif [ ${software} = "shapeit5" ]; then
    set_up_shapeit5
    phase_with_shapeit5
  elif [ ${software} = "eagle2" ]; then
    phase_with_eagle2
  else
    echo "Software ${software} is not valid"
  fi
else
  print_update "Warning: ${out} already exists! Skipping." | tee /dev/stderr
fi

make_tabix ${out}



