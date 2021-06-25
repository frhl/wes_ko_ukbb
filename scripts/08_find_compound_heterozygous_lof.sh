#!/bin/bash

#$ -cwd
# -N tmp_lof
# -o tmp_lof.log
# -e tmp_lof.errors.log
# -q short.qc
#$ -P lindgren.prjc
# -pe shmem 1
# -t 1-22

CHR=${SGE_TASK_ID}

# IN_FILES
VEP_FILE=/well/lindgren/UKBIOBANK/flassen/projects/KO/IMPUTATION/derived/vep/ukbb_mfi_vep_all_HIGH_maf_info_hwe_wb.txt # note - done on WB data
IN_ROOT=/well/lindgren/UKBIOBANK/flassen/projects/KO/IMPUTATION/derived/phased/ensembl_gene_ranges/
IN_FILE=phased_ukb_plof_regions_maf_info_500k_1206v_chr${CHR}.vcf.gz
FAM_FILE=/well/lindgren/UKBIOBANK/stefania/RelGroups/2020_03_10/QC.fam
FAM_FILE_WB=/well/lindgren/UKBIOBANK/stefania/RelGroups/2020_03_10/QCWB.fam

# OUT_FILES
OUT_ROOT_ALL=/well/lindgren/UKBIOBANK/flassen/projects/KO/IMPUTATION/derived/lof/all/compound_heterozygous/
OUT_ROOT_WB=/well/lindgren/UKBIOBANK/flassen/projects/KO/IMPUTATION/derived/lof/wb/compound_heterozygous/
OUT_FILE=ukb_mfi_compound_heterozygous_chr${CHR}.txt
OUT_FILE_MIN=ukb_mfi_compound_heterozygous_simple_chr${CHR}.txt
OUT_FILE_WB=ukb_mfi_wb_compound_heterozygous_chr${CHR}.txt
OUT_FILE_WB_MIN=ukb_mfi_wb_compound_heterozygous_simple_chr${CHR}.txt
TMP_VEP=vep_compound_hetz${CHR}.tmp

module load R

###################
## find variants ##
###################

## subset VEP (chr, snp, gene, INFO)
#cat $VEP_FILE | awk '$1==1' | head -n 20 | tr "|" ";" | awk -F" |;" '{print $1" "$4" "$18" "$13","$14","$15","$16","$17","$18","$19","$20","$22","$26","$27","$28}' | sed 's/CSQ=//' > ${OUT_ROOT_ALL}${TMP_VEP}
cat $VEP_FILE | tr "|" ";" | awk -F" |;" '{print $1" "$4" "$18" "$13","$14","$15","$16","$17","$18","$19","$20","$22","$26","$27","$28}' | sed 's/CSQ=//' > ${OUT_ROOT_ALL}${TMP_VEP}

## run R
Rscript utils/run_get_compound_hetz_ko.R ${IN_ROOT}${IN_FILE} ${OUT_ROOT_ALL}${TMP_VEP} ${OUT_ROOT_ALL}${OUT_FILE}

## subset by white british ancestry
cp ${OUT_ROOT_ALL}${OUT_FILE}  ${OUT_ROOT_WB}${OUT_FILE_WB}.tmp
join -t, -1 1 -2 1  <(cat ${OUT_ROOT_WB}${OUT_FILE_WB}.tmp | sort -k1)\
                    <(cat ${FAM_FILE_WB} | cut -f1 | uniq | sort)\
                    >${OUT_ROOT_WB}${OUT_FILE_WB} 

rm ${OUT_ROOT_WB}${OUT_FILE_WB}.tmp ${OUT_ROOT_ALL}${TMP_VEP}


#########################################################
## reduce them down to simplest informative components ##
#########################################################

## filter files (ALL)
cat ${OUT_ROOT_ALL}${OUT_FILE} | \
  grep COMPOUND_HETEROZYGOUS | \
  cut -d"," -f1,2,3,4,5,6,7,8,9,10 | \
  tr ',' ' ' | \
  awk '{print $1","$3","$5","$6","$2","$9","$8","$4}' >\
  ${OUT_ROOT_ALL}${OUT_FILE_MIN}

## filter files (WB)
cat ${OUT_ROOT_WB}${OUT_FILE_WB} | \
  grep COMPOUND_HETEROZYGOUS | \
  cut -d"," -f1,2,3,4,5,6,7,8,9,10 | \
  tr ',' ' ' | \
  awk '{print $1","$3","$5","$6","$2","$9","$8","$4}' >\
  ${OUT_ROOT_WB}${OUT_FILE_WB_MIN}






