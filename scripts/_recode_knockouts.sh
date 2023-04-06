#!/usr/bin/env bash

source utils/bash_utils.sh
source utils/qsub_utils.sh

echo "$(pwd)"

readonly rscript=${1?Error: Missing arg2 (in_vcf)}
readonly input_path=${2?Error: Missing arg2 (in_vcf)}
readonly vep_path=${3?Error: Missing arg2 (in_vcf)}
readonly out_prefix=${4?Error: Missing arg2 (in_vcf)}

set_up_rpy
Rscript ${rscript} \
  --input_path ${input_path} \
  --vep_path ${vep_path} \
  --out_prefix ${out_prefix}



