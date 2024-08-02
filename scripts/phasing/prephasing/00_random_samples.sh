#!/usr/bin/env bash
#
# @description create matrix table of random samples to be used for read-backed phasing
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=random_samples
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/random_samples.log
#SBATCH --error=logs/random_samples.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 3
#SBATCH --array=22
#
#$ -N random_samples
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/random_samples.log
#$ -e logs/random_samples.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly curwd=$(pwd)
readonly cluster=$( get_current_cluster)
readonly hail_script="scripts/phasing/prephasing/00_random_samples.py"
readonly spark_dir="data/tmp/spark"

readonly input_dir="data/unphased/wes_union_calls/prefilter/200k"
readonly input_path="${input_dir}/ukb_wes_union_calls_chr21.vcf.gz" 
readonly input_type="vcf"

readonly out_dir="data/prephased/wes_union_calls/intervals"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_random_samples_10k_seed1995_chr${chr}"

readonly random_samples_count=10000
readonly seed=1995

mkdir -p ${out_dir}


SECONDS=0
set_up_hail
set_up_pythonpath_legacy
python3 ${hail_script} \
  --random_samples_count ${random_samples_count} \
  --input_path ${input_path} \
  --input_type ${input_type} \
  --out_prefix ${out_prefix} \
  --seed ${seed}



