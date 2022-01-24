#!/usr/bin/env bash
#
#$ -N _summary_statistics_merge
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_summary_statistics_merge.log
#$ -e logs/_summary_statistics_merge.errors.log
#$ -V

source utils/bash_utils.sh

readonly prefix=${1?Error: Missing arg1 (prefix)}
readonly out=${2?Error: Missing arg2 (out)}

# Create header
file=$(echo ${prefix} | sed -e "s/CHR//g")
zcat ${file}* | tail -n+2 > ${out}

# Create remaining content
for chr in {1..22}; do
   file=$(echo ${prefix} | sed -e "s/CHR/${chr}/g")
   zcat "${file}.txt.gz" | tail -n +2  >> ${out}
   echo "concatenating ${file} to ${out}.."
done
gzip ${out}







