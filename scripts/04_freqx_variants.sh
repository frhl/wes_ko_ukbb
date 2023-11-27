#!/usr/bin/env bash
#
# @description Annotate main MatrixTables with VEP results
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=freqx_variants
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/freqx_variants.log
#SBATCH --error=logs/freqx_variants.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --dependency="afterok:38165677"
#SBATCH --array=1-21

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_dir="data/mt/prefilter/pp90"
readonly in_vcf="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.vcf.bgz"

readonly out_dir="data/mt/prefilter/pp90"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005"

mkdir -p ${out_dir}

module load PLINK/1.9b_6.21-x86_64

plink --vcf ${in_vcf} --freqx --out ${out_prefix}

rm ${out_prefix}.log
rm ${out_prefix}.nosex


