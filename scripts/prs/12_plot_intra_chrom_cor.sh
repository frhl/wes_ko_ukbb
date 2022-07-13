#!/usr/bin/env bash
#
#$ -N plot_intra_chrom_cor
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/plot_intra_chrom_cor.log
#$ -e logs/plot_intra_chrom_cor.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qf
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/prs/12_plot_intra_chrom_cor.R"

readonly validation_dir="data/prs/validation/chrom/old"
readonly out_dir="data/prs/validation"
readonly out_prefix="${out_dir}/intra_chrom_correlation_pgs"

set_up_rpy
Rscript "${rscript}" \
    --out_prefix "${out_prefix}" \
    --in_dir "${validation_dir}" \
    --out_height 10 \
    --out_width 12





