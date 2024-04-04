#!/usr/bin/env bash
#
#$ -N aggr_anno
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/aggr_anno.log
#$ -e logs/aggr_anno.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 12
#$ -q short.qc
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/04_aggr_anno.R"

readonly out_dir="data/knockouts/tables"
readonly out_prefix="${out_dir}/combined_annotations.all_syn"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
 --out_prefix "${out_prefix}" \



