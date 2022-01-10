#!/usr/bin/env bash
#
#
#$ -N extract_saige_results
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/extract_saige_results.log
#$ -e logs/extract_saige_results.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#& -V


source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly rscript="scripts/09_extract_saige_results.R"
readonly out_dir="derived/tables/saige"
readonly out_prefix="${out_dir}/211220_wes200k_saige_merge"

# run r-code
mkdir -p ${out_dir}
set_up_rpy

# run for binary traits
set -x
Rscript "${rscript}"\
  --out_prefix ${out_prefix}\
  --in_phenotypes "binary"
set +x

# run for continious traits
set -x
Rscript "${rscript}"\
  --out_prefix ${out_prefix}\
  --in_phenotypes "cts"
set +x



