#!/usr/bin/env bash
#
#$ -N nbglm_reactome
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/nbglm_reactome.log
#$ -e logs/nbglm_reactome.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/12_nbglm_reactome.R"

readonly in_dir="data/knockouts/tables"
readonly in_file="${in_dir}/combined_annotations_by_sample.nohets.txt.gz"

readonly out_dir="data/knockouts/tables"
readonly out_prefix="${out_dir}/nbglm_reactome"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --out_prefix "${out_prefix}" \
  --consolidated "${in_file}"


