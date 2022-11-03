#!/usr/bin/env bash
#
# @description: extract 200K samples 
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=extract_samples
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/extract_samples.log
#SBATCH --error=logs/extract_samples.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=21
#
#
#$ -N extract_samples
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/extract_samples.log
#$ -e logs/extract_samples.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qe
#$ -t 21
#$ -V


source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/11_extract_samples.py"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

## REMOVE TEST PATHS BELOW

readonly in_dir="data/phased/calls/shapeit5/500k"
readonly in_file="${in_dir}/ukb_phased_calls_500k_chr${chr}.vcf.gz"
readonly in_type="vcf"

readonly out_dir="data/phased/calls/shapeit5/200k_from_500k"
readonly out_prefix="${out_dir}/ukb_phased_calls_200k_from_500k_chr${chr}"
readonly out_type="vcf"

readonly samples="data/unphased/overlap/ukb_calls_wes_samples.txt"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

if [ ! -f "${out_prefix}.vcf.bgz" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --input_path "${in_file}" \
     --input_type "${in_type}" \
     --extract_samples "${samples}" \
     --out_prefix "${out_prefix}" \
     --out_type "${out_type}"
fi

module purge
module load BCFtools/1.12-GCC-10.3.0
make_tabix "${out_prefix}.vcf.bgz" "tbi"



