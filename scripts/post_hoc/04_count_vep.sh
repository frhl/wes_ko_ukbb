#!/usr/bin/env bash
#
#$ -N count_vep
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/count_vep.log
#$ -e logs/count_vep.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/post_hoc/04_count_vep.R"

readonly in_dir="data/mt/vep/worst_csq_by_gene_canonical"
readonly in_file_prefilter="${in_dir}/ukb_eur_wes_union_calls_200k_chrCHR.tsv.gz"
readonly in_file_postfilter="${in_dir}/ukb_eur_wes_union_calls_200k_chrCHR.pp90.tsv.gz"
readonly out_dir="data/post_hoc/results"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_variants.pp90"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --in_file_prefilter ${in_file_prefilter} \
  --in_file_postfilter ${in_file_postfilter} \
  --out_prefix ${out_prefix}









