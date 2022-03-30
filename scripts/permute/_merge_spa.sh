#!/usr/bin/env bash
#
#$ -N _merge_spa
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_merge_spa.log
#$ -e logs/_merge_spa.errors.log
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly prefix=${1?Error: Missing arg1 (prefix)}
readonly actual_tasks=${2?Error: Missing arg2 (n_tasks)}
readonly max_tasks=${3?Error: Missing arg2 (n_tasks)}
readonly out_no_gz="${prefix}.txt"

readonly files="${prefix##*/}$"
echo ">${files}"

for id in $(seq 1 ${actual_tasks}); do
   file="${prefix}_${id}of${max_tasks}.txt"
   if [ -f ${file} ]; then
     echo ${file}
     if [ "${id}" == "1" ]; then
        cat "${file}" | head -n 1  >> "${out_no_gz}"
     fi
     cat "${file}" | tail -n +2  >> "${out_no_gz}"
     rm -f "${file}"
   fi
done


if [ -f ${out_no_gz} ]; then
  gzip "${out_no_gz}"
  echo "Merge complete for ${out_no_gz}"
else 
  >&2 "Error: ${out_no_gz} could not be found."
fi




