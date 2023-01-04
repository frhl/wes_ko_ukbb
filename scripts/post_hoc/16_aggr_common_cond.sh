#!/usr/bin/env bash
#
#$ -N aggr_common_cond
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/aggr_common_cond.log
#$ -e logs/aggr_common_cond.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

source utils/bash_utils.sh

readonly rscript="scripts/post_hoc/16_aggr_common_cond.R"

# folder containing ukb_eur_wes_200k_DM_T1D_pLoF_damaging_missense_cond_ENSG00000213676.markers for example
readonly markers_dir="data/conditional/common/spa_iter"
readonly out_dir="data/post_hoc/results"
readonly out_prefix="${out_dir}/aggr_common_cond"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --out_prefix "${out_prefix}" \
  --markers_dir "${markers_dir}"

