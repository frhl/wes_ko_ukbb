#!/usr/bin/env bash
#
#$ -N aggr_prs
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/aggr_prs.log
#$ -e logs/aggr_prs.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/prs/07_aggr_prs.R"

readonly in_dir="data/prs/scores/auto"

readonly pheno_dir="data/phenotypes"
readonly phenos="${pheno_dir}/curated_phenotypes_header.tsv"

readonly out_dir="data/prs/scores"
readonly out_prefix="${out_dir}/ukb_eur_500k_hm3_binary_pgs"

mkdir -p ${out_dir}

aggregate_pgs()
{ 
  set_up_rpy
  set -x
  Rscript "${rscript}" \
   --path_phenos "${phenos}" \
   --in_dir "${in_dir}" \
   --out_prefix "${out_prefix}" 
  set +x
}

aggregate_pgs




