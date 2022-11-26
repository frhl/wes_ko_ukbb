#!/usr/bin/env bash
#

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/conditional/rare/_merge_variant_tables.R"

readonly in_prefix=${1?Error: Missing arg1 (in_prefix)}
readonly chunks=${2?Error: Missing arg2 (chunks)}
readonly out_prefix=${3?Error: Missing arg3 (out_prefix)}

set_up_rpy
Rscript "${rscript}" \
 --in_prefix ${in_prefix} \
 --chunks ${chunks} \
 --out_prefix ${out_prefix} 


