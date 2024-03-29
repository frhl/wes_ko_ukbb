#!/usr/bin/env bash
#
# @description combine whole exome sequences with variants from genotyping array
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=wes_union_calls_gen
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/wes_union_calls_gen.log
#SBATCH --error=logs/wes_union_calls_gen.errors.log
#SBATCH --partition=long
#SBATCH --cpus-per-task 3
#SBATCH --array=20

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/legacy/02_geno_gen.py"
readonly in_dir="data/unphased/wes/post-qc"

readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} )
readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly in_type="mt"

readonly out_dir="data/unphased/wes_union_calls/legacy_long"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}"
readonly out_type="mt"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}
#if [ ! -f "${out_prefix}.vcf.bgz" ]; then
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy
  set -x
  python3 "${hail_script}" \
     --chrom "${chr}" \
     --input_path "${in_file}" \
     --input_type "${in_type}" \
     --out_prefix "${out_prefix}" \
     --out_type "${out_type}" \
     --min_mac 2 \
     --missing 0.05 \
     --dataset "calls" \
     --exclude_trio_parents \
     --export_parents \
     --liftover
  set +x
#else
#  print_update "file ${out} already exists. Skipping!"
#fi

#module purge
#module load BCFtools/1.12-GCC-10.3.0
#make_tabix "${out_prefix}.vcf.bgz" "tbi"

