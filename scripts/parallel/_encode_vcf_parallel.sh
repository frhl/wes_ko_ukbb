#!/usr/bin/env bash

# encoding VCF for interval of genes as opposed
# to an entire chromosome, which avoids Hail out of memory errors

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/_encode_vcf_parallel.py"
readonly spark_dir="data/tmp/spark"

set -x
readonly input_path=${1?Error: Missing arg1 (input_path)}
readonly input_type=${2?Error: Missing arg2 (input_type)}
readonly in_category=${3?Error: Missing arg3 (variant category)}
readonly gene_intervals=${4?Error: Missing arg4 (gene_intervals)}
readonly out_prefix=${5?Error: Missing arg5 (path prefix for saige output)}
readonly genes_per_chunk=${6?Error: Missing arg5 (path prefix for saige output)}
readonly total_chunks=${7?Error: Missing arg5 (path prefix for saige output)}
readonly chr=${7?Error: Missing arg6 (chr)}

readonly chunk_id=$( get_array_task_id )
readonly out_prefix_chunk="${out_prefix}.${chunk_id}of${total_chunks}"
readonly out_prefix_checkpoint="${out_prefix_chunk}_checkpoint.mt/"

readonly idx_start="$(( ((${chunk_id}-1) * ${genes_per_chunk})+1 ))"
readonly idx_end="$(( ${chunk_id} * ${genes_per_chunk}))"
readonly idx_quit=$((${idx_end}+1))
readonly genes=$(cat ${gene_intervals} | sed -n "${idx_start},${idx_end}p;${idx_quit}q" | cut -f2 | tr "\n" ",")
echo ${genes}

evaluate_knockouts() {
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
      --input_path ${input_path} \
      --input_type ${input_type} \
      --csqs_category ${in_category} \
      --out_prefix ${out_prefix_chunk} \
      --subset_gene ${genes}
 rm -rf "${out_prefix_checkpoint}"
}
    

if [ ! -f "${out_prefix_chunk}.tsv.gz" ]; then
  evaluate_knockouts
else
  >&2 echo "${out_prefix_chunk} already exists!"
fi

