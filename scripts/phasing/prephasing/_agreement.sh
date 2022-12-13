#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/vcf_utils.sh

readonly input_path=${1?Error: Missing arg1 (input_path)}
readonly out_prefix=${2?Error: Missing arg3 (out_prefix)}
readonly n_samples=${3?Error: Missing arg3 (n_samples)}
readonly wes_variants=${4?Error: Missing arg4 (wes_variants)}
readonly pp_cutoff=${5?Error: Missing arg5 (pp_cutoff)}
readonly rscript=${6?Error: Missing arg6 (rscript)}


set_up_rpy
Rscript ${rscript} \
  --input_path "${input_path}" \
  --n_samples ${n_samples} \
  --seed 52 \
  --pp_cutoff ${pp_cutoff} \
  --sites "${wes_variants}" \
  --out_prefix "${out_prefix}"


