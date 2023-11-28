#!/usr/bin/env bash

#SBATCH -A lindgren.prj
#SBATCH -J phenotypes_tte
#SBATCH -D /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/phenotypes_tte.log
#SBATCH --error=logs/phenotypes_tte.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1
#
#
#$ -N phenotypes_tte
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/phenotypes_tte.log
#$ -e logs/phenotypes_tte.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly rscript="scripts/00_phenotypes_tte.R"

readonly out_dir="data/phenotypes"
readonly out_path="${out_dir}/tte_matrix_176k.txt.gz"

readonly header_dir="data/phenotypes"
readonly path_header="${header_dir}/dec22_phenotypes_binary_200k_header.tsv"

readonly samples_dir="data/phenotypes/samples"
readonly samples="${samples_dir}/ukb_wes_ko.qc.nfe.samples"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --samples ${samples} \
  --path_header ${path_header} \
  --out_path ${out_path}



