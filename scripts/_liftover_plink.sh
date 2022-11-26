#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly in_prefix=${1?Error: Missing arg1 (in_prefix)}
readonly out_prefix=${2?Error: Missing arg2 (out_prefix)}
readonly in_type=${3?Error: Missing arg3 (in_type)}
readonly out_type=${4?Error: Missing arg4 (out_type)}

readonly chr="${SLURM_ARRAY_TASK_ID}"
readonly in_prefix_chr="${in_prefix}_chr${chr}"
readonly out_prefix_chr="${out_prefix}_chr${chr}"

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/_liftover_plink.py"

if [ ! -f "${out_prefix_chr}.bed" ]; then
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
    --in_prefix ${in_prefix_chr} \
    --out_prefix ${out_prefix_chr} \
    --in_type ${in_type} \
    --out_type ${out_type} \
    --liftover \
    --only_valid_contigs
else
  echo "${out_prefix_chr}.bed already exists. Skipping.."
fi




