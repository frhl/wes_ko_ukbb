#!/usr/bin/env bash
#
# @description combine whole exome sequences with variants from genotyping array
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=wes_union_calls
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/wes_union_calls.log
#SBATCH --error=logs/wes_union_calls.errors.log
#SBATCH --partition=long
#SBATCH --cpus-per-task 5
#SBATCH --array=21
#
#$ -N wes_union_calls
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/wes_union_calls.log
#$ -e logs/wes_union_calls.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qe
#$ -t 21
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/vcf_utils.sh

module load BCFtools/1.12-GCC-10.3.0

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/03_wes_union_calls.py"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_wes_dir="data/unphased/wes/prefilter"
readonly in_wes_file="${in_wes_dir}/ukb_wes_prefilter_200k_chr${chr}.vcf"

readonly in_calls_dir="data/unphased/calls/prefilter"
readonly in_calls_file="${in_calls_dir}/ukb_prefilter_calls_200k_chr${chr}.vcf"

readonly out_dir="data/unphased/wes_union_calls/bcftools"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}"
readonly out_file="${out_prefix}.vcf.gz"

# keep track of wes and calls variants
readonly wes_variants="${out_prefix}_wes.txt"
readonly calls_variants="${out_prefix}_calls.txt"
# save copy of calls file with variants filtered out
readonly tmp_file="${out_prefix}_tmp.vcf.gz"

get_variants() {
  local _vcf=${1}
  local _out=${2}
  bcf
}





mkdir -p ${out_dir}

if [ -f "${out_prefix}.vcf.bgz" ]; then
  




  bcftools merge ${in_wes_file} ${in_calls_file} -oZ -o ${out_file}
  make_tabix "${out_file}" "tbi"
fi


