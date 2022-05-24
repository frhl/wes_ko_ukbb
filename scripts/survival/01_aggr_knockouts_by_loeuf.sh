#!/usr/bin/env bash
#
#$ -N aggr_knockouts_by_loeuf
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/aggr_knockouts_by_loeuf.log
#$ -e logs/aggr_knockouts_by_loeuf.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1
#$ -tc 1
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/survival/01_aggr_knockouts_by_loeuf.R"

readonly in_pattern="pLoF_damaging_missense_all.tsv.gz"
readonly in_dir="data/knockouts/alt"
readonly out_dir="data/survival/knockouts/eur_no_fin/pLoF_damaging_missense"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --in_pattern ${in_pattern} \
  --in_dir ${in_dir} \
  --out_dir ${out_dir}




