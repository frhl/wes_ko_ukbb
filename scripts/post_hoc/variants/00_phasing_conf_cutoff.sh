#!/usr/bin/env bash
#
#$ -N phasing_conf_cutoff
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/phasing_conf_cutoff.log
#$ -e logs/phasing_conf_cutoff.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qc@@short.hga
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/00_phasing_conf_cutoff.R"

readonly in_dir="data/mt/phase_conf"
readonly pattern="ukb_wes_union_calls_200k_chr"

readonly out_dir="data/knockouts/tables"
readonly out_prefix="${out_dir}/shapeit5_phasing_conf_cutoff_table"
readonly chr="${SGE_TASK_ID}"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
 --in_dir ${in_dir} \
 --regex ${pattern} \
 --out_prefix "${out_prefix}"




