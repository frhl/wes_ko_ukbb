#!/usr/bin/env bash
#
# This script will convert all phased vcf files into bed files and subset to high impact variants.
#
# Author: Frederik Lassen (2021-06-25)
#
#$ -N phase_common
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/make_bed.log
#$ -e logs/make_bed.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qe
#$ -t 22
#$ -V

set -o errexit
set -o nounset
module purge 
chr=${SGE_TASK_ID}

# IN FILES
readonly RAW_ROOT="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/phased"
readonly RAW_FILE="/ukb_wes_200k_phased_chr${chr}.1of1.vcf.gz"

readonly FAM_FILE_WB="/well/lindgren/UKBIOBANK/stefania/RelGroups/2020_10_05/QCWB.txt"
#readonly FAM_FILE_WB="/well/lindgren/UKBIOBANK/stefania/RelGroups/2020_03_10/QCWB.fam"

readonly VEP_FILE="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/derived/vep/ukb_wes_200k_vep_HIGH.txt"

# OUT FILES
readonly OUT_ROOT="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/derived/bed/"
readonly OUT_FILE="/ukb_wes_200k_phased_chr${chr}"
readonly TMP_FILE="/ukb_wes_200k_${chr}.tmp"

## setup paths to sofware
plink=/apps/well/plink/2.00a-20170724/plink2

# select VEP variants and keep in temp file
cat ${VEP_FILE} | cut -f3 > ${OUT_ROOT}${TMP_FILE}

## generate bed file with all variants
$plink \
 --vcf ${RAW_ROOT}${RAW_FILE}  \
 --keep ${FAM_FILE_WB} \
 --extract ${OUT_ROOT}${TMP_FILE} \
 --max-maf 0.02 \
 --geno 0.05 \
 --make-bed \
 --out ${OUT_ROOT}${OUT_FILE}

## calculate hwe stats
$plink \
 --bfile ${RAW_ROOT}${RAW_FILE} \
 --keep ${FAM_FILE_WB} \
 --geno 0.05 \
 --hardy \
 --out ${OUT_ROOT}${OUT_FILE}

rm ${OUT_ROOT}${TMP_FILE}

## generate .bed files in which HWE P-values < 10e-7 are excluded
#cat ${OUT_ROOT_BED_WB}${OUT_BED_FILE}.hardy | awk '$10 < 1e-7 {print $2}' > ${OUT_ROOT_BED_WB}${OUT_HWE_EXL_FILE}
#cat ${OUT_ROOT_BED_WB}${OUT_BED_FILE}.hardy | awk '$10 > 1e-7 {print $2}' > ${OUT_ROOT_BED_WB}${OUT_HWE_ICL_FILE}

#plink \
# --bgen ${OUT_ROOT_BGEN}${OUT_BGEN_FILE} 'ref-first' \
# --sample ${OUT_ROOT_BGEN}${OUT_SAMPLE_FILE} \
# --exclude ${OUT_ROOT_BED_WB}${OUT_HWE_EXL_FILE} \
# --keep ${FAM_FILE_WB} \
# --geno 0.05 \
# --make-bed \
# --out ${OUT_ROOT_BED_WB}${OUT_BED_HWE_FILE}

## calculate hwe stats
#plink \
# --bfile ${OUT_ROOT_BED_WB}${OUT_BED_HWE_FILE} \
# --keep ${FAM_FILE_WB} \
# --geno 0.05 \
# --hardy \
# --out ${OUT_ROOT_BED_WB}${OUT_BED_HWE_FILE}


