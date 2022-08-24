#!/usr/bin/env bash
#
#$ -N aggr_tables
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/aggr_tables.log
#$ -e logs/aggr_tables.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hga
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/01_aggr_tables.R"

readonly out_dir="data/knockouts/tables"
readonly out_prefix="${out_dir}/pLoF_damaging_missense_full"

set_up_rpy
set -x
Rscript "${rscript}" \
 --out_prefix "${out_prefix}"
set +x




