#!/usr/bin/env bash
#
# scripts called by simulate_phenotypes.sh

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/bash_utils.sh

readonly hail_script="scripts/simulation/02_simulate_phenotype.py"
readonly spark_dir="data/tmp/spark_dir"

readonly in_prefix=${1?Error: Missing arg1 (phenotype)}
readonly in_type=${2?Error: Missing arg2 (in_vcf)}
readonly h2=${3?Error: Missing arg3 ()}
readonly b=${4?Error: Missing arg3 ()}
readonly pi=${5?Error: Missing arg3 ()}
readonly K=${6?Error: Missing arg3 ()}
readonly seed=${7?Error: Missing arg3 ()}
readonly out_prefix=${8?Error: Missing arg3 ()}

readonly array_idx=$( get_array_task_id )
readonly out_prefix_jid="${out_prefix}_${array_idx}"
readonly seed_jid=$(( ${array_idx} * ${seed}))

mkdir -p ${spark_dir}


if [ ! -f "${out_prefix_jid}_cols.tsv.gz" ]; then
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --in_prefix "${in_prefix}"\
     --in_type "${in_type}" \
     --h2 ${h2} \
     --b ${b} \
     --pi ${pi} \
     --K ${K} \
     --seed ${seed_jid} \
     --out_prefix "${out_prefix_jid}" \
     && print_update "Finished simulating phenotypes for ${in_prefix}" ${SECONDS} \
     || raise_error "Simulating phenotypes for ${in_prefix} failed"
else
  >&2 echo "${out_prefix_jid} already exists!"
fi

