#!/usr/bin/env bash


source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )

readonly gene_table=${1?Error: Missing arg1 (Gene table intervals)}
readonly padding=${2?Error: Missing arg2 (padding to add upstream/downstream of gens)}
readonly phenotype=${3?Error: Missing arg3 (phenotype)}
readonly out_prefix=${4?Error: Missing arg4 (out prefix)}

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/common/04_make_intervals.py"

set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --padding ${padding} \
   --gene_table ${gene_table} \
   --out_prefix ${out_prefix} \
   && print_update "Finished creating intervals for ${out_prefix}" ${SECONDS} \
   || raise_error "Creating intervals for ${out_prefix} failed!"

