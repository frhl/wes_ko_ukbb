#!/usr/bin/env bash
#
# @description phase variants from UKBB whole exome sequencing using SHAPEIT4.2
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=get_singletons
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/get_singletons.log
#SBATCH --error=logs/get_singletons.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

source utils/bash_utils.sh

readonly rscript="scripts/phasing/experimental/01_get_singletons.R"

readonly chr="${SLURM_ARRAY_TASK_ID}"
readonly in_dir="data/unphased/post-qc"
readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr21.mt"
readonly out_dir="data/reads/urv"
readonly out_path="${out_dir}/post_qc_informative_snps"

mkdir -p ${out_dir}





