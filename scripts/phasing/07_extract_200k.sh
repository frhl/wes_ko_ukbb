#!/usr/bin/env bash
#
# @description: extract 200K samples 
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=extract_200k
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/extract_200k.log
#SBATCH --error=logs/extract_200k.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=21-22

#
#$ -N extract_200k
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/extract_200k.log
#$ -e logs/extract_200k.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qe
#$ -t 21-22
#$ -V


source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/07_extract_200k.py"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly samples=""

readonly in_dir="data/phased/calls/shapeit5/by_maf"
readonly in_file="${in_dir}/ukb_prefilter_calls_500k_chr${chr}.mt"
readonly in_type="vcf"

readonly out_dir="data/phased/calls/shapeit5/by_maf"
readonly out_prefix="${out_dir}/ukb_phased_calls_200k_from_500k_chr${chr}"
readonly out_type="vcf"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}
if [ ! -f "${out_prefix}.vcf.bgz" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --input_path "${in_file}" \
     --input_type "${in_type}" \
     --out_prefix "${out_prefix}" \
     --out_type "${out_type}"
fi

module purge
module load BCFtools/1.12-GCC-10.3.0
make_tabix "${out_prefix}.vcf.bgz" "tbi"



