#!/usr/bin/env bash
#
#$ -N geneset_enrichment
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/geneset_enrichment.log
#$ -e logs/geneset_enrichment.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 3
#$ -q short.qc
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/08_geneset_enrichment.R"

readonly out_dir="data/knockouts/tables"
readonly out_prefix="${out_dir}/pLoF_damaging_missense_geneset_enricment"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
 --out_prefix "${out_prefix_variable}" \



