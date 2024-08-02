#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/_recode_knockouts_merge.R"

readonly regex_prefix=${1?Error: Missing arg1 (input_path)}
readonly total_chunks=${2?Error: Missing arg4 (gene_intervals)}
readonly outfile=${3?Error: Missing arg5 (path prefix for saige output)}

readonly directory="$( dirname ${regex_prefix} )"
readonly filename="$( basename ${regex_prefix} )"

readonly suffix=".txt.gz"
readonly expt_files=${total_chunks}
readonly found_files=$( ls ${directory} | grep ${filename} | grep -E "\\.(txt)|(tsv)\\.gz" | wc -l)

if [ ${expt_files} = ${found_files} ]; then
  set_up_rpy
  Rscript ${rscript} \
    --in_dir ${directory} \
    --regex ${filename} \
    --n_expt ${expt_files} \
    --suffix ${suffix} \
    --out ${outfile}
else
  >&2 echo "Error: Expected ${expt_files} genes but found ${found_files} (regex: ${regex_prefix})."
fi

