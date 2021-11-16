#!/usr/bin/env bash
#
#$ -N combine_qc_tables
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/combine_qc_tables.log
#$ -e logs/combine_qc_tables.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q short.qc@@short.hge

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly spark_dir="data/tmp/spark"
readonly out_dir="derived/tables/worst_csq_for_variant_canonical"

# hail script
readonly hail_script="scripts/11_combine_qc_tables.py"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_unphased"

# run hail
mkdir -p ${out_dir}
set_up_hail
set_up_pythonpath_legacy  
python3 "${hail_script}" \
     --out_prefix ${out_prefix}

print_update "Finished running HAIL for chr${chr}" "${SECONDS}"





