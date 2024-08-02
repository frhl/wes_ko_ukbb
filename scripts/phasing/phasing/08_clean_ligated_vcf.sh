#!/usr/bin/env bash
#
# @description combine chunks in a phase-aware manner into full chromosomes.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=clean_ligated
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/clean_ligated.log
#SBATCH --error=logs/clean_ligated.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=1-21

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/vcf_utils.sh

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_dir="data/phased/wes_union_calls/200k/shapeit5/ligated"
readonly in_prefix="${in_dir}/ukb_wes_union_calls_200k_chr${chr}"
readonly in_vcf="${in_prefix}.vcf.bgz"

readonly out_dir="data/phased/wes_union_calls/200k/shapeit5/clean_ligated"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}"
readonly out_vcf="${out_prefix}.vcf.bgz"

mkdir -p ${out_dir}

module load BCFtools
bcftools view -c 1 -c 1:nonmajor ${in_vcf} -Ou | bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%ALT' -Oz -o ${out_vcf}
bcftools index ${out}





