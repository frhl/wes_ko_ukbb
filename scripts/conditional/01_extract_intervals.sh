#!/usr/bin/env bash
#
#
#$ -N extract_intervals
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/extract_intervals.log
#$ -e logs/extract_intervals.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc


source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly rscript="scripts/conditional/01_extract_intervals.R"
readonly in_dir="derived/tables/saige/"
readonly in_prefix="211111"
readonly out_dir="derived/tables/gene_intervals"
readonly out_prefix="${out_dir}/211111_wes200k_saige_wes"

# run r-code
mkdir -p ${out_dir}
set_up_rpy
Rscript "${rscript}" \
  --in_dir ${in_dir} \
  --in_prefix ${in_prefix} \
  --out_prefix ${out_prefix}




