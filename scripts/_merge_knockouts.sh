#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly regex_prefix=${1?Error: Missing arg1 (input_path)}
readonly intervals=${2?Error: Missing arg4 (gene_intervals)}
readonly outfile=${3?Error: Missing arg5 (path prefix for saige output)}

readonly directory="$( dirname ${regex_prefix} )"
readonly filename="$( basename ${regex_prefix} )"

readonly expt_files=$( cat ${intervals} | wc -l )
readonly found_files=$( ls ${directory} | grep ${filename} | wc -l)

if [ ${expt_files} = ${found_files} ]; then
  echo "input code"
fi







