#!/usr/bin/env bash
#
#$ -N gene_phenotype_case_control_chets
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/gene_phenotype_case_control_chets.log
#$ -e logs/gene_phenotype_case_control_chets.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly rscript="scripts/post_hoc/26_gene_phenotype_case_control_chets.R"

readonly co_occurence_dir="data/knockouts/alt/pp90/co_occurence3"
readonly co_occurence_file="${co_occurence_dir}/co_occurence_by_phenotype_chrCHR.txt.gz"


readonly out_dir="data/knockouts/alt/pp90/co_occurencec3"
readonly out_prefix="${out_dir}/co_occurence_collapsed_pLoF_damaging_missense"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --co_occurence_file ${co_occurence_file} \
  --out_prefix ${out_prefix}


