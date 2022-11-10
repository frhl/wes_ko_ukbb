#!/usr/bin/env bash
#
# @description filter WES quality-controlled data.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=extract_quick
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/extract_quick.log
#SBATCH --error=logs/extract_quick.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=21

# --dependency="afterok:8444324"
#
#$ -N extract_quick
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o data/help/logs/extract_quick.log
#$ -e data/help/logs/extract_quick.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -t 21
#$ -V

set -o errexit
set -o nounset

module purge
source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/help/tmp/spark"
readonly hail_script="scripts/phasing/_read_with_hail.py"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_dir="data/phased/wes_union_calls/200k/shapeit5/phase_common"
readonly in_file="${in_dir}/ukb_wes_union_calls_200k_chr${chr}_phase_common.vcf.gz"
readonly in_type="vcf"

readonly out_dir="data/phased/wes_union_calls/200k/shapeit5/phase_common"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}_phase_common_quick"
readonly out_type="vcf"



set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --input_path ${in_file} \
   --input_type ${in_type} \
   --out_prefix "${out_prefix}" \
   --out_type "${out_type}" \
   --min_maf 0.001





