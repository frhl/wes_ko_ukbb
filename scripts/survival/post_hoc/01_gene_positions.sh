#!/usr/bin/env bash

#SBATCH --account=lindgren.prj
#SBATCH --job-name=gene_positions
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/gene_positions.log
#SBATCH --error=logs/gene_positions.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly rscript="scripts/survival/04_gene_positions.R"

readonly in_dir="data/survival/tables/"
readonly in_file="${in_dir}/230609_forest_table.all.txt.gz"

readonly genes="data/genes/220310_ensgid_grch38_pos.tsv.gz"

readonly out_dir="data/survival/gene_positions"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_PHENO_pLoF_damaging_missense"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --coordinates ${genes} \
  --in_file ${in_file} \
  --out_prefix ${out_prefix}





