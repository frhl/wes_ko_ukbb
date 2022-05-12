#!/usr/bin/env bash
##
#$ -N _merge_markers
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/merge_markers.log
#$ -e logs/merge_markers.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 9
#$ -q short.qc

module load BCFtools/1.12-GCC-10.3.0

source utils/vcf_utils.sh
source utils/hail_utils.sh

readonly arg_mt=${1?Error: Missing arg1 (in dir)}
readonly arg_marker_table=${2?Error: Missing arg2 (in dir)}
readonly arg_out_prefix=${3?Error: Missing arg3 (in dir)}

readonly chr="${SGE_TASK_ID}"
readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/common/04_merge_markers.py"

sub_chr () {
  echo ${1} | sed -e "s/CHR/${chr}/g"
}

readonly mt=$(sub_chr ${arg_mt})
readonly marker_table=$(sub_chr ${arg_marker_table})
readonly out_prefix=$(sub_chr ${arg_out_prefix})


merge_markers() {
  SECONDS=0
  set -x
  python3 "${hail_script}" \
     --input_mt "${mt}" \
     --chrom "chr${chr}" \
     --input_markers "${marker_table}" \
     --out_prefix "${out_prefix}"
  set +x
  make_tabix "${out_prefix}.vcf.bgz" "csi"
  print_update "Hail finished writing VCFs in ${SECONDS}"
}

# run analysis
set_up_hail
set_up_pythonpath_legacy
merge_markers


