#!/usr/bin/env bash
#
#$ -N geneset_overlap
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/geneset_overlap.log
#$ -e logs/geneset_overlap.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 6
#$ -q short.qc
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/02_geneset_overlap.R"

readonly out_dir="data/knockouts/tables"
readonly out_prefix="${out_dir}/geneset_overlap"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
 --out_prefix "${out_prefix}" \



