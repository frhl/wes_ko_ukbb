#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly rscript="scripts/conditional/rare/_make_variant_tables.R"

readonly in_vcf=${1?Error: Missing arg1 (in_vcf)}
readonly vcf_lines=${2?Error: Missing arg2 (vcf_lines)}
readonly lines_per_chunk=${3?Error: Missing arg3 (lines_per_chunk)}
readonly total_chunks=${4?Error: Missing arg4 (total_chunks)}
readonly pheno_file=${5?Error: Missing arg5 (pheno_file)}
readonly pheno_list_csv=${6?Error: Missing arg6 (pheno_list_csv)}
readonly covar_path=${7?Error: Missing arg7 (covar_path)}
readonly out_prefix=${8?Error: Missing arg8 (out_prefix)}

readonly cluster=$( get_current_cluster)
readonly chunk=$( get_array_task_id )
readonly out_prefix_chunk="${out_prefix}.${chunk}of${total_chunks}"

readonly hash_file="${out_prefix_chunk}.hash.txt.gz" 
readonly ac_file="${out_prefix_chunk}.AC.txt.gz" 


if [ ! -f "${hash_file}" ]; then
  set_up_rpy
  Rscript "${rscript}" \
    --in_vcf "${in_vcf}" \
    --vcf_lines "${vcf_lines}" \
    --lines_per_chunk "${lines_per_chunk}" \
    --chunk "${chunk}" \
    --total_chunks "${total_chunks}" \
    --pheno_file "${pheno_file}" \
    --phenotypes "${pheno_list_csv}" \
    --covariates "${covar_path}" \
    --out_prefix "${out_prefix_chunk}"
else
  >&2 echo "${hash_file} already exists. Skipping.."
fi 


