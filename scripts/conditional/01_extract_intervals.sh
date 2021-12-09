#!/usr/bin/env bash
#
# Extract genetic coordinates for significant genes in primary analysis.
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

readonly in_dir="derived/tables/saige"
readonly pheno_dir="data/phenotypes"
readonly out_dir="data/conditional/common/extract_intervals"

readonly phenotypes_binary="${pheno_dir}/UKBB_WES200k_binary_phenotypes_header.txt"
readonly phenotypes_cts="${pheno_dir}/UKBB_WES200k_cts_phenotypes_header.txt"

readonly rscript="scripts/conditional/01_extract_intervals.R"

readonly out_prefix_binary="${out_dir}/211111_wes200k_saige_binary_wes"
readonly out_prefix_cts="${out_dir}/211111_wes200k_saige_cts_wes"

mkdir -p ${out_dir}
set_up_rpy

# binary traits
Rscript "${rscript}" \
  --in_dir "${in_dir}" \
  --in_prefix "211111" \
  --in_phenotypes "${phenotypes_binary}" \
  --out_prefix "${out_prefix_binary}"

# cts traits
Rscript "${rscript}" \
  --in_dir "${in_dir}" \
  --in_prefix "211111" \
  --in_phenotypes "${phenotypes_cts}" \
  --out_prefix "${out_prefix_cts}"


