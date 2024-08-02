#!/usr/bin/env bash
#
# @description filter WES quality-controlled data.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=split_wes
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/split_wes.log
#SBATCH --error=logs/split_wes.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 2
#SBATCH --array=1-22

#
#$ -N split_wes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o data/help/logs/split_wes.log
#$ -e data/help/logs/split_wes.errors.log
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
readonly hail_script="scripts/phasing/prefilter/split_parents.py"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_dir="data/unphased/wes/prefilter/200k"
readonly in_file="${in_dir}/ukb_prefilter_wes_200k_chr${chr}.mt"
readonly in_type="mt"

readonly out_dir="data/unphased/wes/prefilter/200k"
readonly out_prefix="${out_dir}/ukb_split_wes_200k_chr${chr}"
readonly out_prefix_no_parents="${out_prefix}_no_parents"
readonly out_prefix_parents="${out_prefix}_parents"
readonly out_type="vcf"

readonly entry_fields_to_drop="GQ,DP,AD,PL"

# samples overlapping exomes and genotypes
readonly parents_dir="data/unphased/overlap"
readonly parents_path="${parents_dir}/ukb_calls_wes_samples_parents.txt"

mkdir -p ${out_dir}
mkdir -p ${spark_dir}

if [ ! -f "${out_prefix_no_parents}.vcf.bgz" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --input_path ${in_file} \
     --input_type ${in_type} \
     --parents_path ${parents_path} \
     --out_prefix "${out_prefix}" \
     --out_type "${out_type}"
fi

if [ ! -f "${out_prefix_no_parents}.vcf.bgz.tbi" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_prefix_no_parents}.vcf.bgz" "tbi"
fi

if [ -f "${out_prefix_parents}.vcf.bgz" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_prefix_parents}.vcf.bgz" "tbi"
fi




