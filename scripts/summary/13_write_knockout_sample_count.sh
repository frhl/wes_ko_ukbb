#!/usr/bin/env bash
#
#$ -N write_knockout_sample_count
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/write_knockout_sample_count.log
#$ -e logs/write_knockout_sample_count.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qa
#$ -V

source utils/bash_utils.sh

# directories
readonly rscript="scripts/summary/13_write_knockout_sample_count.R"

readonly out_prefix="derived/tables/knockouts/ukb_wes200k"

# run hail
set_up_rpy
Rscript "${rscript}" \
     --out_prefix ${out_prefix}






