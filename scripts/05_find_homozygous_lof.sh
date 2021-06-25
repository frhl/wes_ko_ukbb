#!/bin/bash

# finds homozygous LoF variants in the data. Generates a
# result for all samples and white british samples.

#$ -cwd
#$ -P lindgren.prjc
# -N find_homo_lof
# -o find_homo_lof.log
# -e find_homo_lof.errors.log
# -q short.qc
# -pe shmem 5
# -t 1-22

CHR=${SGE_TASK_ID}

# load modules required
module load python/2.7.11


# I/O (WB)
VEP_FILE_WB=/well/lindgren/UKBIOBANK/flassen/projects/KO/IMPUTATION/derived/vep/ukbb_mfi_vep_all_HIGH_maf_info_hwe_wb.txt
IN_ROOT_WB=/well/lindgren/UKBIOBANK/flassen/projects/KO/IMPUTATION/derived/plink/wb/
IN_FILE_WB=ukb_vep_hwe_chr${CHR}.bed
OUT_ROOT_WB=/well/lindgren/UKBIOBANK/flassen/projects/KO/IMPUTATION/derived/lof/wb/homozygous/
OUT_FILE_WB=ukb_mfi_homo_plof_wb_chr${CHR}.txt
OUT_FILE_VEP_WB=ukb_mfi_homo_plof_wb_mrg_chr${CHR}.txt

python utils/get_find_homo_plof.py ${IN_ROOT_WB}${IN_FILE_WB} ${OUT_ROOT_WB}${OUT_FILE_WB}

# merge with VEP consequence (note, that this also contains multi allelic variants).
join -1 5 -2 4 <(cat ${OUT_ROOT_WB}${OUT_FILE_WB} | awk '!/^#/' | sort -k5)  <(cat ${VEP_FILE_WB} | sort -k4) | \
 awk '{print $2" "$3" "$4" "$5" "$1" "$9" "$14" "$18}' | \
 awk -F"|" '{print $1" "$4" "$5" "$6" "$8" "$16" "$15" "$14}' \
 > ${OUT_ROOT_WB}${OUT_FILE_VEP_WB}

