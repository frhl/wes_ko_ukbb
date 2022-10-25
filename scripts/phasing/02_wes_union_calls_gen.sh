#!/usr/bin/env bash
#
# @description combine whole exome sequences with variants from genotyping array
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=wes_union_calls_gen
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/wes_union_calls_gen.log
#SBATCH --error=logs/wes_union_calls_gen.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 5
#SBATCH --array=1-19

#
#$ -N wes_union_calls
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/wes_union_calls.log
#$ -e logs/wes_union_calls.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qe
#$ -t 1-19
#$ -V


source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/02_wes_union_calls_gen.py"
readonly in_dir="data/unphased/wes/prefilter"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )
readonly in_file="${in_dir}/ukb_eur_wes_prefilter_200k_chr${chr}.vcf.bgz"
readonly in_type="vcf"

readonly out_dir="data/unphased/wes_union_calls/new"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}"
readonly out_type="vcf"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}
if [ ! -f "${out_prefix}.vcf.bgz" ]; then
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --chrom "${chr}" \
     --input_path "${in_file}" \
     --input_type "${in_type}" \
     --out_prefix "${out_prefix}" \
     --out_type "${out_type}" \
     --missing 0.05 \
     --dataset "calls" \
     --min_mac 1 \
     --exclude_trio_parents \
     --export_parents \
     --checkpoint \
     --liftover
else
  print_update "file ${out} already exists. Skipping!"
fi

module purge
module load BCFtools/1.12-GCC-10.3.0
make_tabix "${out_prefix}.vcf.bgz" "tbi"



