#!/usr/bin/env bash

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/_merge_plink.py"
readonly spark_dir="data/tmp/spark_dir"

readonly in_prefix=${1?Error: Missing arg2 (in_vcf)}
readonly in_type=${2?Error: Missing arg2 (in_vcf)}
readonly out_prefix=${3?Error: Missing arg2 (in_vcf)}
readonly out_type=${4?Error: Missing arg2 (in_vcf)}

mkdir -p ${out_dir}


if [ ! -f "${out_prefix}.bed" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --in_prefix "${in_prefix}" \
     --in_type "${in_type}" \
     --out_prefix "${out_prefix}" \
     --out_type "${out_type}"
else
  echo "${out_prefix}.bed already exists. Skipping.."
fi


