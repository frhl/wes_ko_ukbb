#!/usr/bin/env bash
#
#$ -N _prs_test
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_prs.log
#$ -e logs/_prs.errors.log
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly r_script=${1?Error: Missing arg1 (r_script)}
readonly pred=${2?Error: Missing arg3 (prediction file)}
readonly ldsc=${3?Error: Missing arg2 (ld_matrix)}
readonly ld_dir=${4?Error: Missing arg2 (ld_matrix)}
readonly method=${5?Error: Missing arg2 (method)}
readonly impute=${6?Error: Missing arg2 (impute)}
readonly prefix=${7?Error: Missing arg8 (prefix)}

readonly chr="${SGE_TASK_ID}"
readonly pred_chr=$(echo ${pred} | sed -e "s/CHR/${chr}/g")
readonly out_prefix_chr=$(echo ${prefix} | sed -e "s/CHR/${chr}/g")
readonly tmp_bfile="${out_prefix_chr}.bfile"

export OPENBLAS_NUM_THREADS=1 # avoid two levels of parallelization

if [ ! -f "${out_prefix_chr}.txt.gz" ]; then
  set_up_ldpred2
  duration=SECONDS
  set -x
  Rscript "${r_script}" \
      --chrom "chr${chr}" \
      --pred "${pred_chr}" \
      --ldsc "${ldsc}" \
      --ld_dir "${ld_dir}" \
      --impute "${impute}" \
      --method "${method}" \
      --tmp_bfile "${tmp_bfile}" \
      --out_prefix "${out_prefix_chr}"
  set +x
  log_runtime $duration
else
  echo "Note: ${out_prefix_chr} already exists. Skipping.."
fi

