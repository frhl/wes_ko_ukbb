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
readonly covariates=${6?Error: Missing arg6 (covariates)}
readonly min_cases=${7?Error: Missing arg6 (covariates)}
readonly prefix=${8?Error: Missing arg8 (out_prefix)}

readonly chr=$( get_array_task_id )
readonly input_path_chr=$(echo ${input_path} | sed -e "s/CHR/${chr}/g")
readonly out_prefix_chr=$(echo ${prefix} | sed -e "s/CHR/${chr}/g")

readonly spark_dir="data/tmp/spark"

echo "nslots: ${NSLOTS}"

set_up_hail 0.2.97
module load OpenBLAS/0.3.1-GCC-7.3.0-2.30
export LD_PRELOAD=/apps/eb/skylake/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib/libopenblas.so

if [ ! -f "${out_prefix_chr}.txt.gz" ]; then
  set -x
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --chrom "${chr}" \
     --input_path "${input_path_chr}" \
     --input_type "${input_type}" \
     --phenotypes "${pheno_file}" \
     --min_cases "${min_cases}" \
     --response "${phenotype}" \
     --covariates "${covariates}" \
     --out_prefix "${out_prefix_chr}" \
     --adjust_maf_by_case_control 
  set +x
else
  echo "Note: ${out_prefix_chr} already exists!"
fi

