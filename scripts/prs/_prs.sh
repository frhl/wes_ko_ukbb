#!/usr/bin/env bash
#
#$ -N _prs
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_prs.log
#$ -e logs/_prs.errors.log
#$ -V

source utils/bash_utils.sh

readonly r_script=${1?Error: Missing arg1 (r_script)}
readonly ld_matrix=${2?Error: Missing arg2 (ld_matrix)}
readonly sumstat=${3?Error: Missing arg2 (sumstat)}
readonly pred=${4?Error: Missing arg3 (prediction file)}
readonly prefix=${5?Error: Missing arg8 (prefix)}

readonly chr="${SGE_TASK_ID}"
readonly ld_matrix_chr=$(echo ${ld_matrix} | sed -e "s/CHR/${chr}/g")
readonly sumstat_chr=$(echo ${sumstat} | sed -e "s/CHR/${chr}/g")
readonly pred_chr=$(echo ${pred} | sed -e "s/CHR/${chr}/g")
readonly out_prefix_chr=$(echo ${prefix} | sed -e "s/CHR/${chr}/g")

set_up_rpy
#if [ ! -f "${out_prefix_chr}.txt.gz" ]; then
set -x
Rscript "${r_script}" \
    --chrom "${chr}" \
    --path_bed_pred "${pred_chr}" \
    --path_ld_matrix "${ld_matrix_chr}" \
    --path_sumstat "${sumstat_chr}" \
    --out_prefix "${out_prefix_chr}"
set +x
#else
#  echo "Note: ${out_prefix_chr} already exists!"
#fi

