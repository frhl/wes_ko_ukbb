#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly hail_script="scripts/06_alleles_by_gene.py"
readonly spark_dir="data/tmp/spark"

readonly input_path=${1?Error: Missing arg1 (input_path)}
readonly input_type=${2?Error: Missing arg2 (input_type)}
readonly in_category=${3?Error: Missing arg3 (variant category)}
readonly out_prefix=${4?Error: Missing arg6 (path prefix for saige output)}

readonly task_id=$( get_array_task_id )
readonly chr=${task_id}

readonly input_path_chr=$(echo ${input_path} | sed -e "s/CHR/${chr}/g")
readonly out_prefix_chr=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")
readonly out="${out_prefix_chr}_alleles.txt.gz"

if [ ! -f "${out}" ]; then
  set +u
  set_up_hail
  set -u
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
      --input_path ${input_path_chr} \
      --input_type ${input_type} \
      --csqs_category ${in_category} \
      --out_prefix ${out_prefix_chr}
else
  >&2 echo "${out} already exists. Skipping"
fi







