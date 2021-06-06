#!/bin/bash

# Takes the resulting VEP files and create two files with AF < 2% 
# that contains HIGH and MODERATE impact variants respectively.

#$ -cwd
# -o combine_vep.log
# -e combine_vep.errors.log
# -q short.qc
#$ -P lindgren.prjc
	# -pe shmem 1

# I/O
IN_ROOT=/well/lindgren/UKBIOBANK/flassen/projects/KO/WES200K/derived/vep/output/
IN_FILES=ukb_wes_200k_inliers_split_filtered_chr*
OUT_ROOT=/well/lindgren/UKBIOBANK/flassen/projects/KO/WES200K/derived/vep/
OUT_FILE1=ukb_wes_200k_inliers_maf_high.txt
OUT_FILE2=ukb_wes_200k_inliers_maf_moderate.txt

# output all high impact files
cat ${IN_ROOT}${IN_FILES} | grep HIGH | awk '$6 < 0.02' > ${OUT_ROOT}${OUT_FILE1}

# print out some stats
echo "HC lines: `cat ${OUT_ROOT}${OUT_FILE1} | grep HC | wc -l` "
echo "LC lines: `cat ${OUT_ROOT}${OUT_FILE1} | grep LC | wc -l` "

# output all high moderate
cat ${IN_ROOT}${IN_FILES} | grep MODERATE | awk '$6 < 0.02' > ${OUT_ROOT}${OUT_FILE2}



