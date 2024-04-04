#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=symlink_tables
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/symlink_tables.log
#SBATCH --error=logs/symlink_tables.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1
#SBATCH --array=1-22

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_dir="$(pwd)/data/conditional/rare/combined/chunks/2024"
readonly in_ac="${in_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense_AC.txt.gz"
readonly in_hash="${in_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense_hash.txt.gz"

readonly out_dir="data/conditional/rare/combined/mt"
readonly out_ac="${out_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense_AC.txt.gz"
readonly out_hash="${out_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense_hash.txt.gz"

if [[ -f "${in_ac}" ]]; then
  echo "Symlinking ${in_ac} -> ${out_ac}"
  ln -sf ${in_ac} ${out_ac}
  ln -sf ${in_hash} ${out_hash}
else
  >&2 echo "${in_ac} does not exist. Skipping!"
fi





