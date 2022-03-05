#!/usr/bin/env bash
#
#$ -N phenotypes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/phenotypes.log
#$ -e logs/phenotypes.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qf
#$ -V

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly hail_script="scripts/00_phenotypes.py"
readonly spark_dir="data/tmp/spark"

readonly in_dir="data/phenotypes"
readonly in_bin="${in_dir}/curated_phenotypes_binary.tsv"
readonly in_cts="${in_dir}/curated_phenotypes_cts.tsv"

readonly out_dir="data/phenotypes"
readonly out_bin="${out_dir}/filtered_phenotypes_binary"
readonly out_cts="${out_dir}/filtered_phenotypes_cts"

readonly final_sample_list="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list"

set_up_hail
set_up_pythonpath_legacy  

python3 "${hail_script}" \
     --input_path "${in_bin}" \
     --extract_samples "${final_sample_list}" \
     --export_header \
     --count_case_control \
     --out_prefix "${out_bin}"

python3 "${hail_script}" \
     --input_path "${in_cts}" \
     --extract_samples "${final_sample_list}" \
     --export_header \
     --out_prefix "${out_cts}"

# create seperate header for Primary Care data
cat "${in_bin}_"
cat "${in_bin}_header.tsv" | grep care > "${out_bin}_PC_header.tsv"
cat "${in_bin}_header.tsv" | grep -v care > "${out_bin}_notPC_header.tsv"

# create seperate header for non residuals (i.e. 
# the signal that's left after conditioning
cat "${in_cts}" | grep residual > "${out_cts}_residual.tsv"


