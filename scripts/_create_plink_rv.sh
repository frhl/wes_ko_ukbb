#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly out_prefix=${1?Error: Missing arg1 (in_vcf)}
readonly samples_keep=${2?Error: Missing arg2 (in_csi)}
readonly rare_markers_per_chrom=${3?Error: Missing arg3}
readonly low_freq_markers_per_chrom=${4?Error: Missing arg4}

# setup out paths
readonly chr="${SLURM_ARRAY_TASK_ID}"
readonly out_prefix_counts="${out_prefix}_counts_chr${chr}"
readonly out="${out_prefix}_chr${chr}"

# path to hard coded genotypes
readonly geno_dir="/well/lindgren/UKBIOBANK/DATA/CALLS"
readonly calls_bed="${geno_dir}/ukb_cal_chr${chr}_v2.bed"
readonly calls_bim="${geno_dir}/ukb_snp_chr${chr}_v2.bim"
readonly calls_fam="/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"


if [ ! -f "${out}.bed" ]; then
  module load PLINK/2.00a2.3_x86_64
  #1. calculate allele counts for each marker in the large plink file with hard called genotypes
  plink2 \
    --bed ${calls_bed} \
    --bim ${calls_bim} \
    --fam ${calls_fam} \
    --keep ${samples_keep} \
    --freq counts \
    --out ${out_prefix_counts}

  #2 randomly extract IDs for markers falling in the two MAC categories
  cat <(tail -n +2 ${out_prefix_counts}.acount | awk '(($6-$5) < 20 && ($6-$5) >= 10) || ($5 < 20 && $5 >= 10) {print $2}' | shuf -n ${rare_markers_per_chrom}) \
    <(tail -n +2 ${out_prefix_counts}.acount | awk ' $5 >= 20 && ($6-$5)>= 20 {print $2}' | shuf -n ${low_freq_markers_per_chrom}) > ${out_prefix_counts}.markerid.list

  # 3. extract markers from the large plink file
  plink2 \
    --bed ${calls_bed} \
    --bim ${calls_bim} \
    --fam ${calls_fam} \
    --keep ${samples_keep} \
    --extract ${out_prefix_counts}.markerid.list \
    --make-bed \
    --out ${out}
  else
  echo "${out} already exists. Skipping.."
fi

