#!/usr/bin/env bash
#
#$ -N calc_intra_chrom_cor
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/calc_intra_chrom_cor.log
#$ -e logs/calc_intra_chrom_cor.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/prs/11_calc_intra_chrom_cor.R"

readonly pgs_dir="data/prs/scores"
readonly out_dir="data/prs/validation/chrom"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
    --pgs_dir "${pgs_dir}" \
    --out_dir "${out_dir}"


