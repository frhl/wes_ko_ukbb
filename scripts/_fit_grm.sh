#!/usr/bin/env bash

source utils/bash_utils.sh

readonly createSparseGRM="utils/saige/createSparseGRM.R"

readonly plink_file=${1?Error: Missing arg1 (plink_file)}
readonly out_prefix=${2?Error: Missing arg2 (out_prefix)}
readonly threads="42"

SECONDS=0
set_up_RSAIGE
set -x
Rscript "${createSparseGRM}" \
    --plinkFile="${plink_file}" \
    --nThreads="${threads}" \
    --outputPrefix="${out_prefix}" \
    --numRandomMarkerforSparseKin=2000 \
    --relatednessCutoff=0.05 \
    && print_update "Finished fitting GRM" ${SECONDS} \
    || raise_error "Fitting GRM failed"
set +x


