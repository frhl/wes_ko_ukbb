#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly rscript="scripts/simulation/_merge_phenotype.R"

readonly in_prefix=${1?Error: Missing arg1 (phenotype)}
readonly pheno_file=${2?Error: Missing arg1 (phenotype)}
readonly covariates=${3?Error: Missing arg1 (phenotype)}
readonly out_file=${4?Error: Missing arg3 ()}

SECONDS=0
set_up_rpy
Rscript "${rscript}" \
   --input_path "${in_prefix}" \
   --real_phenotype_path "${pheno_file}" \
   --covars_keep "${covariates}" \
   --output_path "${out_file}" \
   && print_update "Finished simulating phenotypes for ${in_prefix}" ${SECONDS} \
   || raise_error "Simulating phenotypes for ${in_prefix} failed"

