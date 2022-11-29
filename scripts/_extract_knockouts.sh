#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/_extract_knockouts.py"
readonly spark_dir="data/tmp/spark"

readonly input_path=${1?Error: Missing arg1 (input_path)}
readonly input_type=${2?Error: Missing arg2 (input_type)}
readonly in_category=${3?Error: Missing arg3 (variant category)}
readonly gene_intervals=${4?Error: Missing arg4 (gene_intervals)}
readonly out_prefix=${5?Error: Missing arg5 (path prefix for saige output)}
readonly chr=${6?Error: Missing arg6 (chr)}

readonly task_id=$( get_array_task_id )
readonly gene=$( cat ${gene_intervals} | sed "${task_id}q;d" | cut -f2)
readonly out_prefix_gene="${out_prefix}_${gene}"

evaluate_knockouts() {
  echo "Evaluating ${out_prefix_gene}.vcf.bgz.."
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
      --input_path ${input_path} \
      --input_type ${input_type} \
      --csqs_category ${in_category} \
      --subset_gene ${gene} \
      --out_prefix ${out_prefix_gene} \
      && print_update "Finished evaluation knockouts for chr${chr}-${gene}" ${SECONDS} \
      || raise_error "Evaluating knockouts for chr${chr}-${gene} failed"
}


if [ ! -f "${out_prefix_gene}.tsv.gz" ]; then
  evaluate_knockouts
else
  >&2 echo "${out_prefix_gene} already exists!"
fi

