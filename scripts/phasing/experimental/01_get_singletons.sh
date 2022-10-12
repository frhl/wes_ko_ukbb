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
readonly in_dir="data/knockouts/alt"
readonly in_file="${in_dir}/ukb_eur_wes_200k_chrCHR_maf0to5e-2_pLoF_damaging_missense_all.tsv.gz"
readonly out_dir="data/reads/singletons"
readonly out_path="${out_dir}/samples_with_damaging_missense_singletons.txt"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
   --in_prefix ${in_file} \
   --out_path ${out_path}




