#!/usr/bin/env bash
#
# Find individuals that are homozygous of for selected variants.
#
# Author: Frederik Lassen (2021-06-25)
#
#$ -N homozgygous
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/find_homozygous.log
#$ -e logs/find_homozygous.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qe
#$ -t 22
#$ -V

set -o errexit
set -o nounset
module purge 
chr=${SGE_TASK_ID}

# load modules required
module load python/2.7.11

# infiles
readonly RAW_ROOT="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/derived/bed"
readonly RAW_FILE="/ukb_wes_200k_phased_chr{chr}.bed"
readonly VEP_FILE="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/derived/vep/ukb_wes_200k_vep_HIGH.txt"

# outfiles 
readonly OUT_ROOT="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/derived/alleles/homozygous"
readonly OUT_TMP="/ukb_wes_200k_phased_homozygous_python_chr${chr}.tmp"
readonly OUT_FILE="/ukb_wes_200k_phased_homozygous_chr${chr}.txt"

# find samples/variants that are homozygoys
python utils/get_find_homo_plof.py ${IN_ROOT_WB}${IN_FILE_WB} ${OUT_ROOT_WB}${OUT_FILE_WB}

# merge with VEP consequence (note, that this also contains multi allelic variants).
#join -1 5 -2 4 <(cat ${OUT_ROOT_WB}${OUT_FILE_WB} | awk '!/^#/' | sort -k5)  <(cat ${VEP_FILE_WB} | sort -k4) | \
# awk '{print $2" "$3" "$4" "$5" "$1" "$9" "$14" "$18}' | \
# awk -F"|" '{print $1" "$4" "$5" "$6" "$8" "$16" "$15" "$14}' \
# > ${OUT_ROOT_WB}${OUT_FILE_VEP_WB}


