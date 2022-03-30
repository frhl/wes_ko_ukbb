#!/usr/bin/env bash
#
#$ -N post_hoc
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/post_hoc.log
#$ -e logs/post_hoc.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hga
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/01_aggr_tables.R"

readonly out_dir="data/knockout/combined"
readonly out_prefix="${out_dir}/pLoF_damaging_missense_full"

set_up_rpy
set -x
Rscript "${rscript}" \
 --out_prefix "${out_prefix}"
set +x




