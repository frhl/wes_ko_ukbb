#!/usr/bin/env bash
#
# @description merge phased vcf with unphased parents for later switch error calculation.
# @note - takes about ~ 24h with 1 a core
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=filter_to_trios
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/filter_to_trios.log
#SBATCH --error=logs/filter_to_trios.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=21
#
#$ -N filter_to_trios
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/filter_to_trios.log
#$ -e logs/filter_to_trios.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -t 21
#$ -V

set -o errexit
set -o nounset

readonly hail_script="scripts/phasing/phasing/05_filter_to_trios.py"

source utils/vcf_utils.sh
source utils/hail_utils.sh
source utils/qsub_utils.sh

#readonly hail_script="scripts/phasing/05_parents.py"
readonly spark_dir="data/tmp/spark"
readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

# fam file for calculating switch errors
readonly pedigree_dir="/well/lindgren/UKBIOBANK/nbaya/resources"
readonly pedigree="${pedigree_dir}/ukb11867_pedigree.fam"
# for shapeit5 phasing
readonly phased_dir="data/phased/wes_union_calls/200k/shapeit5/parents_improved"
readonly phased_path="${phased_dir}/ukb_wes_union_calls_200k_shapeit5_parents_chr${chr}.vcf.gz"
readonly phased_type="vcf"
readonly out_dir="data/phased/wes_union_calls/200k/shapeit5/trios"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_trios_chr${chr}"
# out paths and types
readonly out_vcf="${out_prefix}.vcf.gz"
readonly out_type="vcf"

mkdir -p ${out_dir}

echo "test.."

set_up_hail
set_up_pythonpath_legacy
python3 ${hail_script} \
  --in_file ${phased_path} \
  --in_type ${phased_type}\
  --out_prefix ${out_prefix} \
  --out_type ${out_type}



