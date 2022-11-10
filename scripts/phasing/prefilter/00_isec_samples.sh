#!/usr/bin/env bash
#
# @description find samples overlapping wes and calls
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=isec_samples
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/isec_samples.log
#SBATCH --error=logs/isec_samples.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#
#$ -N isec_samples
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/isec_samples.log
#$ -e logs/isec_samples.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qe
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/prefilter/00_isec_samples.py"
readonly chr="21"

# output samples
readonly out_dir="data/unphased/overlap"
readonly out_prefix="${out_dir}/ukb_calls_wes_samples"
# get samples in whole exome seq data
readonly in_wes_dir="data/unphased/wes/post-qc"
readonly in_wes_file="${in_wes_dir}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly in_wes_type="mt"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --chrom "${chr}" \
   --in_path "${in_wes_file}" \
   --in_type "${in_wes_type}" \
   --out_prefix "${out_prefix}" \
   --remove_withdrawn \
   --dataset "calls"



