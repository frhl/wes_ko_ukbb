#!/usr/bin/env bash
##
#$ -N _merge_markers
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_merge_markers.log
#$ -e logs/_merge_markers.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc

module load BCFtools/1.12-GCC-10.3.0

source utils/vcf_utils.sh
source utils/hail_utils.sh

readonly mt=${1?Error: Missing arg1 (in dir)}
readonly marker_table=${2?Error: Missing arg2 (in dir)}
readonly out_prefix=${3?Error: Missing arg3 (in dir)}

readonly chr="${SGE_TASK_ID}"
readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/04_merge_markers.py"

merge_marker() {
  SECONDS=0
  python3 "${hail_script}" \
     --mt ${mt} \
     --chrom ${chr} \
     --input_markers ${marker_table} \
     --out_prefix ${out_prefix} \
  print_update "Hail finished writing VCFs in ${SECONDS}"
}

set_up_hail
set_up_pythonpath_legacy
filter_mergers


