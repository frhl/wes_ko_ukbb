#!/usr/bin/env bash
#
#$ -N _prs_test
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_prs_by_chrom.log
#$ -e logs/_prs_by_chrom.errors.log
#$ -V

source utils/bash_utils.sh

readonly r_script=${1?Error: Missing arg1 (r_script)}
readonly gwas=${2?Error: Missing arg2 (sumstat)}
readonly pred=${3?Error: Missing arg3 (prediction file)}
readonly ld_bed=${4?Error: Missing arg2 (ld_matrix)}
readonly ld_dir=${5?Error: Missing arg2 (ld_matrix)}
readonly method=${6?Error: Missing arg2 (ld_matrix)}
readonly trait=${7?Error: Missing arg2 (ld_matrix)}
readonly prefix=${8?Error: Missing arg8 (prefix)}

readonly chr="${SGE_TASK_ID}"
readonly pred_chr=$(echo ${pred} | sed -e "s/CHR/${chr}/g")
readonly out_prefix_chr=$(echo ${prefix} | sed -e "s/CHR/${chr}/g")

set_up_rpy
if [ ! -f "${out_prefix_chr}.txt.gz" ]; then
  set -x
  Rscript "${r_script}" \
      --chrom "chr${chr}" \
      --gwas "${gwas}" \
      --pred "${pred_chr}" \
      --ld_bed "${ld_bed}" \
      --ld_dir "${ld_dir}" \
      --method "${method}" \
      --trait "${trait}" \
      --out_prefix "${out_prefix_chr}"
  set +x
else
  echo "Note: ${out_prefix_chr} already exists. Skipping.."
fi

