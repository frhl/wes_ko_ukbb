#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly out_prefix=${1?Error: Missing arg1 (out_prefix)}
readonly samples_keep=${2?Error: Missing arg2 (samples_keep)}

# setup out paths
readonly chr="${SLURM_ARRAY_TASK_ID}"
readonly out_prefix_chr="${out_prefix}_chr${chr}"

# path to hard coded genotypes
readonly geno_dir="/well/lindgren/UKBIOBANK/DATA/CALLS"
readonly calls_bed="${geno_dir}/ukb_cal_chr${chr}_v2.bed"
readonly calls_bim="${geno_dir}/ukb_snp_chr${chr}_v2.bim"
readonly calls_fam="/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"


if [ ! -f "${out_prefix_chr}.bed" ]; then
  module load PLINK/2.00a2.3_x86_64
  plink2 \
    --bed ${calls_bed} \
    --bim ${calls_bim} \
    --fam ${calls_fam} \
    --keep ${samples_keep} \
    --indep-pairwise 50 5 0.05 \
    --make-bed \
    --out ${out_prefix_chr}
else 
  echo "${out_prefix_chr} already exists. Skipping.."
fi


