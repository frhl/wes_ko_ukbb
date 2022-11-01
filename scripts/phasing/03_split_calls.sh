#!/usr/bin/env bash
#
# @description generate files of genotyped calls
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=prefilter_calls
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/prefilter_calls.log
#SBATCH --error=logs/prefilter_calls.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 3
#SBATCH --array=1-22
#
#
#$ -N prefilter_calls
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prefilter_calls.log
#$ -e logs/prefilter_calls.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 3
#$ -q short.qe
#$ -t 21-22
#$ -V


source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/02_prefilter_calls.py"

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly out_dir="data/unphased/calls/prefilter/new_by_maf"
readonly out_prefix="${out_dir}/ukb_prefilter_calls_200k_chr${chr}"
readonly out_prefix_parents="${out_prefix}_parents"
readonly out_type="vcf"

# samples overlapping exomes and genotypes
readonly samples_dir="data/unphased/overlap"
readonly samples_list="${samples_dir}/ukb_calls_wes_samples.txt"
readonly trio_parents="${samples_dir}/ukb_calls_wes_samples_parents.txt"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

if [ ! -f "${out_prefix}.vcf.bgz" ]; then
  set_up_hail 0.2.97
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --chrom "${chr}" \
     --out_prefix "${out_prefix}" \
     --out_type "${out_type}" \
     --exclude_trio_parents "${trio_parents}" \
     --filter_incorrect_reference \
     --liftover \
     --export_parents \
     --min_maf 0.001 \
     --missing 0.05 \
     --dataset "calls"
fi

if [ ${out_type} == "vcf" ] & [ ! -f "${out_prefix}.vcf.bgz.tbi" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_prefix}.vcf.bgz" "tbi"
fi

if [ -f "${out_prefix_parents}.vcf.bgz" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_prefix_parents}.vcf.bgz" "tbi"
fi





