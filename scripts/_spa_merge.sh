#!/usr/bin/env bash
#
#$ -N _spa_merge
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_spa_merge.log
#$ -e logs/_spa_merge.errors.log
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly prefix=${1?Error: Missing arg1 (prefix)}
readonly out_dir=${2?Error: Missing arg2 (out_dir)}
readonly out=${3?Error: Missing arg3 (file out full name)}
readonly remove_by_chr=${4?Error: Missing arg4 (Remove indiviudal chromosome files)}
readonly out_without_gz=$(echo ${out} | sed -e "s/\\.gz//g" | sed -e "s/_chrCHR//g")

file=$(echo ${prefix} | sed -e "s/CHR/[0-9]+/g" )
files="${file##*/}$"
readonly n=$(ls -l "${out_dir}" | grep -E "${files}" | wc -l)
readonly N=22
# always expecting 22 autosomes
if (( $(echo "$n == $N" | bc -l) )); then
  for chr in {1..22}; do
     file=$(echo ${prefix} | sed -e "s/CHR/${chr}/g")
     if [ "${chr}" == "1" ]; then
        cat "${file}" | head -n 1  >> "${out_without_gz}"
     fi
     cat "${file}" | tail -n +2  >> "${out_without_gz}"
     if [ ${remove_by_chr} == "Y" ]; then
        rm "${file}"
        rm "${file}.index"
     fi 
  done
  echo "Merge completed for ${out_without_gz}."
  gzip "${out_without_gz}"
 else
  >&2 echo "Some chromosomes are missing for ${file} (found ${n} but expected ${N})."
fi





