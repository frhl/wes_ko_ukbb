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
readonly out=${2?Error: Missing arg2 (file out full name)}
readonly remove_by_chr=${3?Error: Missing arg3 (Remove indiviudal chromosome files)}

# set up directories 
readonly out_dir=$( dirname ${out})
readonly out_without_gz=$(echo ${out} | sed -e "s/\\.gz//g" | sed -e "s/_chrCHR//g")

# what files are we looking for
file=$(echo ${prefix} | sed -e "s/CHR/[0-9]+/g" )
files="${file##*/}$"

# count how many are present versus expected
readonly n=$(ls -l "${out_dir}" | grep -E "${files}" | grep -v ".index" | wc -l)
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
        rm -f "${file}"
        rm -f "${file}.index"
     fi 
  done
  echo "Merge completed for ${out_without_gz}."
  gzip "${out_without_gz}"
  rm -f "${out_without_gz}" # remove existing file
 else
  >&2 echo "Some chromosomes are missing for ${file} (found ${n} but expected ${N})."
fi





