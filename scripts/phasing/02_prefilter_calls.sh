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
#SBATCH --array=20-22

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/02_prefilter_calls.py"

readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} )
readonly out_dir="data/unphased/calls/prefilter/by_maf"
readonly out_prefix="${out_dir}/ukb_prefilter_calls_200k_chr${chr}"
readonly out_prefix_parents="${out_prefix}_parents"
readonly out_type="vcf"

# samples overlapping exomes and genotypes
readonly samples_list="data/unphased/overlap/ukb_calls_wes_samples.txt"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

if [ ! -f "${out_prefix}.vcf.bgz" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --chrom "${chr}" \
     --out_prefix "${out_prefix}" \
     --out_type "${out_type}" \
     --extract_samples "${samples_list}" \
     --filter_incorrect_reference \
     --liftover \
     --exclude_trio_parents \
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





