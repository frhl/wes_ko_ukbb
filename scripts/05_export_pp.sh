#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=export_pp
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/export_pp.log
#SBATCH --error=logs/export_pp.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-19
#
#$ -N export_pp
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/export_pp.log
#$ -e logs/export_pp.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 20-22
#$ -V

set -o errexit
set -o nounset

# Note: PP-field is only defiend for MAF < 0.01%

source utils/vcf_utils.sh
source utils/hail_utils.sh
source utils/qsub_utils.sh

readonly hail_script="scripts/05_export_pp.py"
readonly spark_dir="data/tmp/spark"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_dir="data/mt/annotated"
readonly in_path="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.mt"
readonly in_type="mt"

readonly out_dir="data/mt/phase_conf"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}"

mkdir -p ${out_dir}

module purge
set_up_hail
set_up_pythonpath_legacy
python3 ${hail_script} \
  --phased_path ${in_path} \
  --phased_type ${in_type} \
  --out_prefix ${out_prefix}

