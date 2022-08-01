#!/usr/bin/env bash
#
#$ -N _make_variant_tables
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_make_variant_tables.log
#$ -e logs/_make_variant_tables.errors.log
#$ -P lindgren.prjc
#$ -q short.qa
#$ -pe shmem 1
#$ -V

set -o errexit
set -o nounset

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

readonly chunk=${SGE_TASK_ID}
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


