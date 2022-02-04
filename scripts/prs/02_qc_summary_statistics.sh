#!/usr/bin/env bash
#
#$ -N summary_statistics
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/summary_statistics.log
#$ -e logs/summary_statistics.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hga
#$ -t 34
#$ -V

module purge
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly pheno_dir="data/phenotypes"
readonly in_dir="data/prs/sumstat"

readonly r_script="scripts/prs/_summary_statistics.sh"

readonly pheno_file="${pheno_dir}/curated_phenotypes.tsv" 
readonly pheno_list="${pheno_dir}/curated_phenotypes_cts_header.tsv"
readonly phenotype=$( cut -f${SGE_TASK_ID} ${pheno_list} )

readonly out_prefix="${out_dir}/ukb_hapmap_500k_${phenotype}"





