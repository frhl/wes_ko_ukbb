#!/usr/bin/env bash
#
# @description merge phased vcf with unphased parents for later switch error calculation.
# @note - takes about ~ 24h with 1 a core
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=extract_pp_field
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/extract_pp_field.log
#SBATCH --error=logs/extract_pp_field.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=20,22
#SBATCH --dependency="afterok:8613568"
#
#$ -N extract_pp_field
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/extract_pp_field.log
#$ -e logs/extract_pp_field.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -t 21
#$ -V

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/hail_utils.sh
source utils/qsub_utils.sh

readonly hail_script="scripts/phasing/15_extract_pp_field.py"
readonly spark_dir="data/tmp/spark"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly phased_dir="data/phased/wes_scaffold_calls/200k_from_500k/ligated"
readonly phased_path="${phased_dir}/ukb_wes_scaffold_calls_200k_fromm_500k_chr${chr}.vcf.bgz"
readonly phased_type="vcf"

readonly out_dir="data/phased/wes_scaffold_calls/200k_from_500k/ligated"
readonly out_prefix="${out_dir}/ukb_wes_scaffold_calls_200k_from_500k_shapeit5_chr${chr}.quality"

mkdir -p ${out_dir}


# combine parents and children in same vcf
# note: we can't combine data using bcftools beacuse some variants
# share the same position (but not same reference alleles)
module purge
set_up_hail
set_up_pythonpath_legacy
python3 ${hail_script} \
  --phased_path ${phased_path} \
  --phased_type ${phased_type} \
  --out_prefix ${out_prefix} \
  --max_mac "50"

