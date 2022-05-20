#!/usr/bin/env bash
#
#$ -N ldsc_summary
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/ldsc_summary.log
#$ -e logs/ldsc_summary.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qf
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/prs/07_ldsc_summary.R"

readonly ldsc_dir="data/prs/ldsc"
readonly out_dir="data/prs/scores"
readonly out_prefix="${out_dir}/ldsc_summary"

set_up_rpy
Rscript "${rscript}" \
    --in_dir "${ldsc_dir}" \
    --out_prefix "${out_prefix}"


