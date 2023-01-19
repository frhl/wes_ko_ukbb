#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=calc_intra_chrom_cor
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/calc_intra_chrom_cor.log
#SBATCH --error=logs/calc_intra_chrom_cor.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/prs/11_calc_intra_chrom_cor.R"

readonly pgs_dir="data/prs/scores_full"
readonly out_dir="data/prs/validation/230118_chrom"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
    --pgs_dir "${pgs_dir}" \
    --out_dir "${out_dir}"


