#!/usr/bin/env bash
#
#$ -N co_occurence
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/co_occurence.log
#$ -e logs/co_occurence.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -t 1-22
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly rscript="scripts/post_hoc/01_co_occurence.R"
readonly out_dir="data/knockouts/alt/pp90/co_occurence3"

readonly out_prefix="${out_dir}/co_occurence_chr${chr}_other_missense"
readonly annotation="other_missense"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --chrom ${chr} \
  --annotation ${annotation} \
  --out_prefix ${out_prefix}


