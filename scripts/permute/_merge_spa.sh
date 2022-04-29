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

readonly shuffles=${1?Error: Missing arg1 (prefix)}
readonly replicates=${2?Error: Missing arg1 (prefix)}
readonly phenotype=${3?Error: Missing arg1 (prefix)}
readonly out_prefix=${4?Error: Missing arg2 (n_tasks)}
readonly out=${5?Error: Missing arg2 (n_tasks)}

aggregate_saige() {

  local prefix="${out_prefix}_${phenotype}"
  local out_no_gz="${out}"
  local max_tasks=$(( (${shuffles} / ${replicates})  ))
  rm -f "${out_no_gz}.gz"
  if [ ${shuffles} -ge 100 ]; then
    if [ ${max_tasks} -ge 1 ]; then
      for id in $(seq 1 ${max_tasks}); do
        file="${prefix}_${id}.txt.gz"
         if [ -f ${file} ]; then
           echo ${file}
           if [ "${id}" == "1" ]; then
              zcat "${file}" | head -n 1  >> "${out_no_gz}"
           fi
           zcat "${file}" | tail -n +2  >> "${out_no_gz}"
         else
           >&2 echo "File ${file} does not exists (aggregate_saige). Skipping.."
         fi
       done
       gzip "${out_no_gz}"
       rm -f ${out}
       echo "Aggregated ${max_tasks} files to ${out_no_gz}."
    else
      raise_error "Need at least one task to submit merge"
    fi
  else
    raise_error "Invalid amount of shuffles"
  fi
}

aggregate_saige


