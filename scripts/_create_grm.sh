#!/usr/bin/env bash
#
#$ -N _create_grm
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_create_grm.log
#$ -e logs/_create_grm.errors.log
#$ -P lindgren.prjc

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly out_prefix=${1?Error: Missing arg2 (in_vcf)}
readonly out_type=${2?Error: Missing arg2 (in_vcf)}
readonly final_sample_list=${3?Error: Missing arg3 (in_csi)}

readonly chr="${SGE_TASK_ID}"
readonly out="${out_prefix}_chr${chr}"
readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/_create_grm.py"

SECONDS=0
set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
  --chroms ${chr} \
  --out_prefix ${out} \
  --out_type ${out_type} \
  --final_sample_list ${final_sample_list} \
  --use_markers_by_kinship \
  --use_markers_by_mac 200 \
  && print_update "Finished writing samples for relatedness matrix (GRM)" ${SECONDS} \
  || raise_error "Writing samples for GRM failed"


