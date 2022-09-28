#!/usr/bin/env bash
#
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=plot_intra_chrom_cor
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/plot_intra_chrom_cor.log
#SBATCH --error=logs/plot_intra_chrom_cor.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --requeue

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





