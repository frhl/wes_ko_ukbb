#!/usr/bin/env bash
#
# @description phase variants from UKBB whole exome sequencing using SHAPEIT4.2
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=smartphase
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/smartphase.log
#SBATCH --error=logs/smartphase.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=21

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

readonly in_dir="data/mt/annotated"
readonly var_dir="data/reads/singletons"
readonly out_dir="data/reads/phased"

mkdir -p ${out_dir}

readonly chr="${SLURM_ARRAY_TASK_ID}"
readonly vcf_file="${in_dir}/ukb_eur_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly out_file="${out_dir}/ukb_eur_wes_union_calls_200k_chr${chr}_smartphase.tsv"
readonly var_file="${var_dir}/samples_with_damaging_missense_singletons.txt"


# all variants in probabilistic knockouts
# allele counts of variants



