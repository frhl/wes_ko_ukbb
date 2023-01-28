#!/usr/bin/env bash
#
#$ -N co_occurence_by_phenotype
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/co_occurence_by_phenotype.log
#$ -e logs/co_occurence_by_phenotype.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1-22
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly rscript="scripts/post_hoc/23_co_occurence_by_phenotype.R"

readonly path_header="data/phenotypes/dec22_phenotypes_binary_200k_header.tsv"
readonly path_phenotypes="data/phenotypes/dec22_phenotypes_binary_200k.tsv.gz"

readonly time_to_event_dir="/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources"
readonly time_to_event_samvida="${time_to_event_dir}/eid_time_to_event_matrix.txt"
readonly time_to_event_duncan="${time_to_event_dir}/eid_time_to_event_matrix_duncan_phenotypes.txt"

#readonly out_dir="data/knockouts/alt/pp90/co_occurence_time_to_event"
#readonly out_prefix="${out_dir}/co_occurence_samvida_tte_by_phenotype_chr${chr}"
readonly out_dir="data/knockouts/alt/pp90/co_occurence3"
readonly out_prefix="${out_dir}/co_occurence_by_phenotype_chr${chr}"
readonly annotation="pLoF_damaging_missense"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --chrom ${chr} \
  --annotation ${annotation} \
  --path_phenotypes ${path_phenotypes} \
  --path_header ${path_header} \
  --out_prefix ${out_prefix}


