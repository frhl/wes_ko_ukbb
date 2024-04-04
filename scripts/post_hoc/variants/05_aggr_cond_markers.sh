#!/usr/bin/env bash
#
#$ -N aggr_cond_markers
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/aggr_cond_markers.log
#$ -e logs/aggr_cond_markers.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

source utils/bash_utils.sh

readonly rscript="scripts/post_hoc/05_aggr_cond_markers.R"

# folder containing ukb_eur_wes_200k_DM_T1D_pLoF_damaging_missense_cond_ENSG00000213676.markers for example
readonly markers_common_dir="data/conditional/common/spa_iter"

# folder containing rare markers used, e.g. ukb_eur_wes_200k_chr6_PSOR_combined_pLoF_damaging_missense_locoprs.rare.markers.keep
readonly markers_rare_dir="data/saige/output/binary//step2_combined//min_mac4"

readonly out_dir="data/post_hoc/results"
readonly out_prefix="${out_dir}/aggr_cond_markers"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --out_prefix "${out_prefix}" \
  --markers_common_dir "${markers_common_dir}" \
  --markers_rare_dir "${markers_rare_dir}"

