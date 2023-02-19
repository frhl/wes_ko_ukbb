#!/usr/bin/env bash
#
#$ -N gene_geneset_counts
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/gene_geneset_counts.log
#$ -e logs/gene_geneset_counts.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/24_gene_geneset_counts.R"

readonly out_dir="data/knockouts/tables"
readonly out_prefix="${out_dir}/gene_geneset_counts"

mkdir -p ${out_dir}


get_gene_table () {
  local anotation=${1}
  Rscript "${rscript}" \
    --annotation "${annotation}" \
    --out_prefix "${out_prefix}"
}

set_up_rpy
get_gene_table "pLoF"


