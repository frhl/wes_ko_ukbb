#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly hail_script=${1?Error: Missing arg1 (hail_script)}
readonly input_path=${2?Error: Missing arg2 (input_path)}
readonly input_type=${3?Error: Missing arg2 (input_type)}
readonly pheno_file=${4?Error: Missing arg3 (pheno_file)}
readonly phenotype=${5?Error: Missing arg4 (phenotype)}
readonly covariates=${6?Error: Missing arg5 (covariates)}
readonly prefix=${7?Error: Missing arg7 (out_prefix)}

readonly chr=$( get_array_task_id )
readonly input_path_chr=$(echo ${input_path} | sed -e "s/CHR/${chr}/g")
readonly out_prefix_chr=$(echo ${prefix} | sed -e "s/CHR/${chr}/g")

readonly spark_dir="data/tmp/spark"

set_up_hail_debug
set_up_pythonpath_legacy
export LD_PRELOAD="/well/lindgren/users/mmq446/conda/skylake/envs/hail-v0.2.105/lib/libopenblas.so"

set -x

if [ ! -f "${out_prefix_chr}.txt.gz" ]; then
  python3 "${hail_script}" \
     --chrom "${chr}" \
     --input_path "${input_path_chr}" \
     --input_type "${input_type}" \
     --phenotypes "${pheno_file}" \
     --response "${phenotype}" \
     --covariates "${covariates}" \
     --out_prefix "${out_prefix_chr}"
else
  echo "Note: ${out_prefix_chr}.txt.gz already exists!"
fi

