#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
set_up_hail
set_up_pythonpath_legacy

readonly chr=${1?Error: Missing arg1 (chr)} # Chromosome, e.g. "1" for chrom 1
readonly vcf_to_phase=${2?Error: Missing arg2 (vcf_to_phase)} # Path of VCF to phase
readonly min_interval_unit=${3?Error: Missing arg3 (min_interval_unit)} 
readonly interval_path=${4?Error: Missing arg4 (min_interval_unit)} 
readonly phasing_region_size=${5?Error: Missing arg5 (phasing_region_size)} # Minimum guaranteed size of phasing window in terms of variant count
readonly phasing_region_overlap=${6?Error: Missing arg6 (phasing_region_overlap)} # Minimum overlap between adjacent phasing windows
readonly max_phasing_region_size=${7?Error: Missing arg7 (max_phasing_region_size)} # Maximum size of phasing window allowed, only used at the end of a chromosome. Must be larger than phasing_region_size
readonly out_prefix=${8?Error: Missing arg8 (path prefix for output intermediate VCF)} # Path to output phased VCF
readonly pedigree=${9?Error: Missing arg9 (pedigree)} # Path of VCF to use as haplotype scaffold
readonly software=${10?Error: Missing arg10 (software)} # Path of VCF to use as haplotype scaffold
readonly threads=$(( ${SLURM_CPUS_ON_NODE} - 1))

readonly hail_script="scripts/phasing/04_phase_chunks.py"
readonly interval_flags="--chrom ${chr} --min_interval_unit ${min_interval_unit} --phasing_region_size ${phasing_region_size} --phasing_region_overlap ${phasing_region_overlap} --max_phasing_region_size ${max_phasing_region_size}"
readonly phasing_idx=${SLURM_ARRAY_TASK_ID} # one-based index for which phasing interval to phase

readonly max_phasing_idx=$( python3 ${hail_script} ${interval_flags} --get_max_phasing_idx --interval_path ${interval_path})
readonly out_prefix_w_phasing_idx="${out_prefix}.${phasing_idx}of${max_phasing_idx}"
readonly out="${out_prefix_w_phasing_idx}.vcf.gz"
readonly log="${out_prefix_w_phasing_idx}.log"

phase_with_shapeit() {
  SECONDS=0
  local gmap="/well/lindgren/flassen/software/SHAPEIT4/b38.gmap/chr${chr}.b38.gmap.gz"
  mkdir -p $( dirname ${out} )
  readonly region=$( python3 ${hail_script} ${interval_flags} --get_interval --phasing_idx ${phasing_idx} --interval_path ${interval_path} )
  print_update "Starting SHAPEIT4 phasing for ${region}, out: ${out}"
  set -x
  shapeit4.2 \
    --input ${vcf_to_phase} \
    --map ${gmap} \
    --region "chr${region}" \
    --thread ${threads} \
    --sequencing \
    --log ${log} \
    --output ${out}
  set +x
  duration=${SECONDS}
  print_update "Finished phasing variants for ${region} (chr${chr} ${phasing_idx}/${max_phasing_idx}, out: ${out}" "${duration}"
}

phase_with_eagle2() {
  SECONDS=0
  local eagle="/well/lindgren/flassen/software/eagle/Eagle_v2.4.1/./eagle"
  local gmap="/well/lindgren/flassen/software/eagle/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
  mkdir -p $( dirname ${out} )
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
    phase_with_shapeit  
  elif [ ${software} = "eagle2" ]; then
    phase_with_eagle2
  else
    echo "Software ${software} is not valid"
  fi
else
  print_update "Warning: ${out} already exists! Skipping." | tee /dev/stderr
fi

make_tabix ${out}



