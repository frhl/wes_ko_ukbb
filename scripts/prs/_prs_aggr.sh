#!/usr/bin/env bash

source utils/bash_utils.sh

readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly in_dir=${2?Error: Missing arg2 (out_dir)}
readonly out_dir=${3?Error: Missing arg3 (mrg_dir)}
readonly rscript="scripts/prs/_prs_aggr.R"

set_up_rpy
Rscript ${rscript} \
  --phenotype ${phenotype} \
  --in_dir ${in_dir} \
  --out_dir ${out_dir} \
  --regex ""




