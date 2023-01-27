#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly cluster=$( get_current_cluster )
readonly index=$( get_array_task_id )
readonly chr=$( get_chr ${index} )

readonly pheno_file=${1?Error: Missing arg1 (pheno_file)}
readonly phenotype=${2?Error: Missing arg2 (phenotype)}
readonly annotation=${3?Error: Missing arg3 (phenotype)}
readonly out_prefix=${4?Error: Missing arg4 (out_prefix)}

readonly rscript="scripts/post_hoc/24_co_occurence.R"

set_up_rpy
Rscript ${rscript} \
  --chrom ${chr} \
  --phenotype ${phenotype} \
  --path_phenotypes ${pheno_file} \
  --annotation ${annotation} \
  --out_prefix ${out_prefix}


