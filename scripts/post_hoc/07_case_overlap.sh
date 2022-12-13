#!/usr/bin/env bash
#
#$ -N case_overlap
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/case_overlap.log
#$ -e logs/case_overlap.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 3
#$ -q short.qc
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/post_hoc/07_case_overlap.R"

readonly pheno_dir="data/phenotypes"
readonly phenos="${pheno_dir}/spiros_brava_phenotypes_binary_200k.tsv.gz"
readonly header="${pheno_dir}/spiros_brava_phenotypes_binary_200k_header.tsv"
readonly out_dir="data/post_hoc/results"
readonly out_prefix="${out_dir}/case_control_overlap"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --phenos ${phenos} \
  --header ${phenos} \
  --out_prefix ${out_prefix}





