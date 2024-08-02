#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly chr=${1?Error: Missing arg1 (chr)} # Chromosome, e.g. "1" for chrom 1
readonly vcf_to_phase=${2?Error: Missing arg2 (vcf_to_phase)} # Path of VCF to phase
readonly vcf_to_scaffold=${3?Error: Missing arg3 (vcf_to_scaffold)} # Path of VCF to phase
readonly min_interval_unit=${4?Error: Missing arg4 (min_interval_unit)} 
readonly interval_path=${5?Error: Missing arg5 (interval_path)} 
readonly phasing_region_size=${6?Error: Missing arg6 (phasing_region_size)} # Minimum guaranteed size of phasing window in terms of variant count
readonly phasing_region_overlap=${7?Error: Missing arg7 (phasing_region_overlap)} # Minimum overlap between adjacent phasing windows
readonly max_phasing_region_size=${8?Error: Missing arg8 (max_phasing_region_size)} # Maximum size of phasing window allowed, only used at the end of a chromosome. Must be larger than phasing_region_size
readonly out_prefix=${9?Error: Missing arg9 (out_prefix)} # Path to output phased VCF
readonly software=${10?Error: Missing arg10 (software)} # Path of VCF to use as haplotype scaffold
readonly pbwt_min_mac=${11?Error: Missing arg11 (pbwt_min_mac)} #
readonly pbwt_depth=${12?Error: Missing arg12 (pbwt_depth)} #
readonly pbwt_modulo=${13?Error: Missing arg13 (pbwt_modulo)} #
readonly pbwt_mdr=${14?Error: Missing arg14 (pbwt_ndr)} #
readonly ps_error=${15?Error: Missing arg15 (ps_error)} #
readonly effective_size=${16?Error: Missing arg16 (pop_effective_size)} #
readonly threads=$( get_threads )


readonly hail_script="scripts/phasing/phasing/02_phase_chunks.py"
readonly interval_flags="--chrom ${chr} --min_interval_unit ${min_interval_unit} --phasing_region_size ${phasing_region_size} --phasing_region_overlap ${phasing_region_overlap} --max_phasing_region_size ${max_phasing_region_size}"
readonly phasing_idx=$( get_array_task_id ) # one-based index for which phasing interval to phase

set_up_hail
set_up_pythonpath_legacy
readonly max_phasing_idx=$( python3 ${hail_script} ${interval_flags} --get_max_phasing_idx --interval_path ${interval_path})
readonly out_prefix_w_phasing_idx="${out_prefix}.${phasing_idx}of${max_phasing_idx}"
readonly out="${out_prefix_w_phasing_idx}.vcf.gz"
readonly log="${out_prefix_w_phasing_idx}.log"

mkdir -p $( dirname ${out} )

phase_with_shapeit5() {
  SECONDS=0
  readonly gmap="/well/lindgren/flassen/software/SHAPEIT4/b38.gmap/chr${chr}.b38.gmap.gz"
  readonly region=$( python3 ${hail_script} ${interval_flags} --get_interval --phasing_idx ${phasing_idx} --interval_path ${interval_path} )
  print_update "Starting SHAPEIT5 (rare) phasing for ${region} (pbwt_min_mac >= ${pbwt_min_mac}) out: ${out}"
  ${SHAPEIT_phase_rare} \
    --input-region "chr${region}" \
    --input-plain ${vcf_to_phase} \
    --scaffold ${vcf_to_scaffold} \
    --scaffold-region "chr${region}" \
    --effective-size ${effective_size} \
    --pbwt-mac ${pbwt_min_mac} \
    --pbwt-depth-rare ${pbwt_depth} \
    --pbwt-depth-common ${pbwt_depth} \
    --pbwt-modulo ${pbwt_modulo} \
    --pbwt-mdr ${pbwt_mdr} \
    --map ${gmap} \
    --thread ${threads} \
    --log ${log} \
    --output ${out} \
    && print_update "Finished phasing variants for chr${chr} using SHAPEIT5, out: ${out}" "${SECONDS}" \
    || raise_error "$( print_update "Phasing variants failed for chr${chr}" ${SECONDS} )" 
}

phase_with_shapeit4() {
  SECONDS=0
  readonly gmap="/well/lindgren/flassen/software/SHAPEIT4/b38.gmap/chr${chr}.b38.gmap.gz"
  readonly region=$( python3 ${hail_script} ${interval_flags} --get_interval --phasing_idx ${phasing_idx} --interval_path ${interval_path} )
  print_update "Starting SHAPEIT4 phasing for ${region} (min_mac >= ${pbwt_min_mac}) out: ${out}"
  shapeit4.2 \
    --input ${vcf_to_phase} \
    --map ${gmap} \
    --region "chr${region}" \
    --thread ${threads} \
    --pbwt-mac ${pbwt_min_mac} \
    --pbwt-depth ${pbwt_depth} \
    --log ${log} \
    --output ${out} \
    --sequencing \
    --use-PS ${ps_error} \
    && print_update "Finished phasing variants for chr${chr} using SHAPEIT4, out: ${out}" "${SECONDS}" \
    || raise_error "$( print_update "Phasing variants failed for chr${chr}" ${SECONDS} )" 
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
  local duration=${SECONDS}
  print_update "Finished phasing for ${region} (chr${chr} ${phasing_idx}/${max_phasing_idx}, out: ${out}" "${duration}"
}



echo "Checking if ${out} exists.."
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

#module purge
#module load BCFtools/1.12-GCC-10.3.0
#make_tabix ${out}



