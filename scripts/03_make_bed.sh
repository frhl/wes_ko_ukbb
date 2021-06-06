#!/bin/bash

# genereates .bgen and .bed files
# filters by the following (include if):
# 1) HWE P-value > 1e-7
# 3) genotyping rate > 5%

#$ -cwd
#$ -P lindgren.prjc
#$ -

CHR=${SGE_TASK_ID}

# IN FILES
readonly RAW_ROOT=/well/lindgren/UKBIOBANK/nbaya/resources/ukb_wes_200k_inliers_split_filtered/
readonly RAW_FILE=ukb_wes_200k_inliers_split_filtered_hail_chr${CHR}.vcf.bgz

#readonly SAM_ROOT=/gpfs1/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/
#readonly SAM_FILE=ukb11867_imp_chr1_v3_s487395.sample

readonly FAM_FILE=/well/lindgren/UKBIOBANK/stefania/RelGroups/2020_03_10/QC.fam
readonly FAM_FILE_WB=/well/lindgren/UKBIOBANK/stefania/RelGroups/2020_03_10/QCWB.fam

readonly VEP_FILE=/well/lindgren/UKBIOBANK/flassen/projects/KO/IMPUTATION/derived/vep/ukbb_mfi_vep_all_HIGH_maf_info.txt

# OUT FILES
OUT_ROOT_BGEN=/well/lindgren/UKBIOBANK/flassen/projects/KO/IMPUTATION/derived/plink/bgen/
OUT_ROOT_BED=/well/lindgren/UKBIOBANK/flassen/projects/KO/IMPUTATION/derived/plink/all/
OUT_ROOT_BED_WB=/well/lindgren/UKBIOBANK/flassen/projects/KO/IMPUTATION/derived/plink/wb/
OUT_SAMPLE_FILE=ukb_vep_chr${CHR}.sample
OUT_BGEN_FILE=ukb_vep_chr${CHR}.bgen
OUT_BED_FILE=ukb_vep_chr${CHR}
OUT_BED_HWE_FILE=ukb_vep_hwe_chr${CHR}
OUT_HWE_EXL_FILE=ukb_vep_hwe_chr_excluded_snps_chr${CHR}.txt
OUT_HWE_ICL_FILE=ukb_vep_hwe_chr_included_snps_chr${CHR}.txt
TMP_FILE=rsid_chr${CHR}.tmp

## setup paths to sofware
qctool=/apps/well/qctool/2.0.1/qctool
plink=/apps/well/plink/2.00a-20170724/plink2

## Open VEP variants.
cat ${VEP_FILE} | cut -d" " -f3 > ${OUT_ROOT_BGEN}${TMP_FILE}

## subset SNPs
#$qctool \
# -g ${RAW_ROOT}${RAW_FILE} \
# -s ${SAM_ROOT}${SAM_FILE} \
# -incl-snpids ${OUT_ROOT_BGEN}${TMP_FILE} \
# -ofiletype bgen \
# -og ${OUT_ROOT_BGEN}${OUT_BGEN_FILE} \
# -os ${OUT_ROOT_BGEN}${OUT_SAMPLE_FILE}

# remove tmp files
rm ${OUT_ROOT_BGEN}${TMP_FILE}

###############
# All samples #
###############

## generate .bed files in which SNPs with a genotype rate of 95% (5 % missing) are included.
$plink \
 --bgen ${OUT_ROOT_BGEN}${OUT_BGEN_FILE} 'ref-first' \
 --sample ${OUT_ROOT_BGEN}${OUT_SAMPLE_FILE} \
 --keep ${FAM_FILE} \
 --geno 0.05 \
 --make-bed \
 --out ${OUT_ROOT_BED}${OUT_BED_FILE}

## calculate hwe stats
$plink \
 --bfile ${OUT_ROOT_BED}${OUT_BED_FILE} \
 --geno 0.05 \
 --hardy \
 --out ${OUT_ROOT_BED}${OUT_BED_FILE}

## generate .bed files in which HWE P-values < 10e-7 are excluded
cat ${OUT_ROOT_BED}${OUT_BED_FILE}.hardy | awk '$10 < 1e-7 {print $2}' > ${OUT_ROOT_BED}${OUT_HWE_EXL_FILE}
cat ${OUT_ROOT_BED}${OUT_BED_FILE}.hardy | awk '$10 > 1e-7 {print $2}' > ${OUT_ROOT_BED}${OUT_HWE_ICL_FILE}

$plink \
 --bgen ${OUT_ROOT_BGEN}${OUT_BGEN_FILE} 'ref-first' \
 --sample ${OUT_ROOT_BGEN}${OUT_SAMPLE_FILE} \
 --exclude ${OUT_ROOT_BED}${OUT_HWE_EXL_FILE} \
 --keep ${FAM_FILE} \
 --geno 0.05 \
 --make-bed \
 --out ${OUT_ROOT_BED}${OUT_BED_HWE_FILE}

$plink \
 --bfile ${OUT_ROOT_BED}${OUT_BED_HWE_FILE} \
 --geno 0.05 \
 --hardy \
 --out ${OUT_ROOT_BED}${OUT_BED_HWE_FILE}

#################
# White British #
#################

$plink \
 --bgen ${OUT_ROOT_BGEN}${OUT_BGEN_FILE} 'ref-first' \
 --sample ${OUT_ROOT_BGEN}${OUT_SAMPLE_FILE} \
 --keep ${FAM_FILE_WB} \
 --geno 0.05 \
 --make-bed \
 --out ${OUT_ROOT_BED_WB}${OUT_BED_FILE}

## calculate hwe stats
$plink \
 --bfile ${OUT_ROOT_BED_WB}${OUT_BED_FILE} \
 --keep ${FAM_FILE_WB} \
 --geno 0.05 \
 --hardy \
 --out ${OUT_ROOT_BED_WB}${OUT_BED_FILE}

## generate .bed files in which HWE P-values < 10e-7 are excluded
cat ${OUT_ROOT_BED_WB}${OUT_BED_FILE}.hardy | awk '$10 < 1e-7 {print $2}' > ${OUT_ROOT_BED_WB}${OUT_HWE_EXL_FILE}
cat ${OUT_ROOT_BED_WB}${OUT_BED_FILE}.hardy | awk '$10 > 1e-7 {print $2}' > ${OUT_ROOT_BED_WB}${OUT_HWE_ICL_FILE}

$plink \
 --bgen ${OUT_ROOT_BGEN}${OUT_BGEN_FILE} 'ref-first' \
 --sample ${OUT_ROOT_BGEN}${OUT_SAMPLE_FILE} \
 --exclude ${OUT_ROOT_BED_WB}${OUT_HWE_EXL_FILE} \
 --keep ${FAM_FILE_WB} \
 --geno 0.05 \
 --make-bed \
 --out ${OUT_ROOT_BED_WB}${OUT_BED_HWE_FILE}

## calculate hwe stats
$plink \
 --bfile ${OUT_ROOT_BED_WB}${OUT_BED_HWE_FILE} \
 --keep ${FAM_FILE_WB} \
 --geno 0.05 \
 --hardy \
 --out ${OUT_ROOT_BED_WB}${OUT_BED_HWE_FILE}


