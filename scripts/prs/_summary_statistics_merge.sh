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

file=$(echo ${prefix} | sed -e "s/CHR//g")

#x=$(ls -l data/prs/sumstat | grep 'ukb_hapmap_500k_eur_Alanine_aminotransferase_residual_chr[0-9]\+\.txt\.gz' | wc -l)

readonly files='${file}[0-9]\+\.txt\.gz'
readonly n=$(ls -l ${files} | grep ${files} | wc -l)
readonly N=21

if (( $(echo "$n > $N" | bc -l) )); then
  zcat ${file}* | tail -n+2 > ${out}
  for chr in {1..22}; do
     file=$(echo ${prefix} | sed -e "s/CHR/${chr}/g")
     zcat "${file}.txt.gz" | tail -n +2  >> ${out}
     echo "concatenating ${file} to ${out}.."
  done
  gzip ${out}
else
  >&2 echo "Some chromosomes are missing for ${file}!"
fi





