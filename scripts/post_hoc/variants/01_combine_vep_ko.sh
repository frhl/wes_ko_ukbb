#!/usr/bin/env bash
#
#$ -N combine_vep_ko
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/combine_vep_ko.log
#$ -e logs/combine_vep_ko.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hga
#$ -t 1-22
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/01_combine_vep_ko.R"

readonly out_dir="data/knockouts/tables"
readonly out_prefix="${out_dir}/pLoF_damaging_missense_full"
readonly chr="${SGE_TASK_ID}"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
 --chrom ${chr} \
 --out_prefix "${out_prefix}"




