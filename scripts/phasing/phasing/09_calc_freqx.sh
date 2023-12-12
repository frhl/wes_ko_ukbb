#!/usr/bin/env bash
#
# @description merge phased vcf with unphased calc_freqx for later switch error calculation.
# @note - takes about ~ 24h with 1 a core
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=calc_freqx
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/calc_freqx.log
#SBATCH --error=logs/calc_freqx.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-22

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/hail_utils.sh
source utils/qsub_utils.sh

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly phased_dir="data/phased/wes_union_calls/200k/shapeit5/filter_by_pp"
readonly phased_path="${phased_dir}/ukb_wes_union_calls_200k_chr${chr}.clean.vcf.bgz"
readonly out_dir="data/phased/wes_union_calls/200k/shapeit5/filter_by_pp"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}.clean"

mkdir -p ${out_dir}

module load PLINK/1.9b_6.21-x86_64

plink --vcf ${phased_path} --freqx --out ${out_prefix} && gzip -f ${out_prefix}.frqx
rm ${out_prefix}.log
rm ${out_prefix}.nosex





