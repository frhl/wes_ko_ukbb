#!/usr/bin/env bash
#
#$ -N write_gene_ko
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/write_gene_ko.log
#$ -e logs/write_gene_ko.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 21
#$ -tc 1
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly rscript="scripts/survival/01_write_gene_ko.R"

readonly in_dir="data/knockouts/alt/pp90/combined"
readonly in_file="${in_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense_all.tsv.gz"

readonly out_dir="data/survival/knockouts/pLoF_damaging_missense/chr${chr}"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --in_file ${in_file} \
  --out_prefix ${out_prefix}




