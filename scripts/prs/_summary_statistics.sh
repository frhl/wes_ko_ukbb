#!/usr/bin/env bash
#
#$ -N _summary_statistics
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_summary_statistics.log
#$ -e logs/_summary_statistics.errors.log
#$ -V

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly hail_script=${1?Error: Missing arg1 (hail_script)}
readonly input_path=${2?Error: Missing arg2 (input_path)}
readonly input_type=${3?Error: Missing arg2 (input_type)}
readonly pheno_file=${4?Error: Missing arg3 (pheno_file)}
readonly phenotype=${5?Error: Missing arg4 (phenotype)}
readonly covariates=${6?Error: Missing arg6 (covariates)}
readonly prefix=${7?Error: Missing arg8 (out_prefix)}

readonly chr=${SGE_TASK_ID}
readonly input_path_chr=$(echo ${input_path} | sed -e "s/CHR/${chr}/g")
readonly out_prefix_chr=$(echo ${prefix} | sed -e "s/CHR/${chr}/g")

readonly spark_dir="data/tmp/spark"

set_up_hail
set_up_pythonpath_legacy
module load OpenBLAS/0.3.8-GCC-9.2.0 # required for linear regression
export LD_PRELOAD=/apps/eb/skylake/software/OpenBLAS/0.3.8-GCC-9.2.0/lib/libopenblas.so # required for linear regression
if [ ! -f "${out_prefix_chr}.txt.gz" ]; then
  set -x
  python3 "${hail_script}" \
     --chrom "${chr}" \
     --input_path "${input_path_chr}" \
     --input_type "${input_type}" \
     --phenotypes "${pheno_file}" \
     --response "${phenotype}" \
     --covariates "${covariates}" \
     --out_prefix "${out_prefix_chr}"
  set +x
else
  echo "Note: ${out_prefix_chr} already exists!"
fi

