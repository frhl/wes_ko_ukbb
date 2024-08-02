#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=symlink_knockouts
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/symlink_knockouts.log
#SBATCH --error=logs/symlink_knockouts.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#
#
#$ -N symlink_knockouts
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/symlink_knockouts.log
#$ -e logs/symlink_knockouts.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh

readonly cluster=$( get_current_cluster)

readonly in_fast_dir="$(pwd)/data/knockouts/alt/pp90/fast"
readonly in_collect_dir="$(pwd)/data/knockouts/alt/pp90/collected"
readonly in_extract_dir="$(pwd)/data/knockouts/alt/pp90/extracted"
readonly in_extract_array_dir="$(pwd)/data/knockouts/alt/pp90/extracted_array"
readonly out_dir="$(pwd)/data/knockouts/alt/pp90/combined"
readonly prefix="ukb_eur_wes_200k"

# symlink knockouts where only vcf encoding has been created
ln -fs ${in_fast_dir}/${prefix}* ${out_dir}/.

# symlink the knockouts that have been created on a gene
# by gene basis (and could not be created throughout the 
# previous step due to memory contraints).
ln -fs ${in_extract_dir}/${prefix}* ${out_dir}/.
ln -fs ${in_extract_array_dir}/${prefix}* ${out_dir}/.

# symlink all knockouts that have been collected
# using a single script in hail
ln -fs ${in_collect_dir}/${prefix}* ${out_dir}/.
















