#!/usr/bin/env bash

source utils/bash_utils.sh

readonly prefix=${1?Error: Missing arg1 (prefix)}
readonly out_dir=${2?Error: Missing arg2 (out_dir)}
readonly out=${3?Error: Missing arg3 (out without)}
readonly out_without_gz=$(echo ${out} | sed -e "s/\\.gz//g")

file=$(echo ${prefix} | sed -e "s/CHR//g" )
readonly basename="${file##*/}"
readonly files="${basename}[0-9]+\.txt\.gz"
readonly n=$(ls -l "${out_dir}" | grep -E "${files}" | wc -l)
readonly N=21

if (( $(echo "$n > $N" | bc -l) )); then
  for chr in {1..22}; do
     file=$(echo ${prefix} | sed -e "s/CHR/${chr}/g")
     if [ "${chr}" == "1" ]; then
        zcat "${file}.txt.gz" | head -n 1  >> "${out_without_gz}"
     fi
     zcat "${file}.txt.gz" | tail -n +2  >> "${out_without_gz}"
  done
  gzip "${out_without_gz}"
else
  >&2 echo "Some chromosomes are missing for ${file}"
fi





