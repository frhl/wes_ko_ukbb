#!/usr/bin/env bash
#
# @description evaluate phasing quality stratified by MAF
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=switch
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/switch.log
#SBATCH --error=logs/switch.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=21

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_dir="data/phased/wes_scaffold_calls/200k_from_500k/merge_parents"
readonly in_file="${in_dir}/ukb_wes_scaffold_calls_200k_from_500k_shapeit5_chr${chr}.vcf.gz"

readonly out_dir="data/phased/wes_scaffold_calls/200k_from_500k/merge_parents"
readonly out_file="${out_dir}/ukb_wes_scaffold_calls_200k_from_500k_shapeit5_chr${chr}.switch"

readonly pedigree_dir="/well/lindgren/UKBIOBANK/nbaya/resources"
readonly pedigree="${pedigree_dir}/ukb11867_pedigree.fam"

readonly region="chr${chr}"

mkdir -p ${out_dir}

set_up_shapeit5
${SHAPEIT_switch} \
  --estimation ${in_file} \
  --validation ${in_file} \
  --pedigre ${pedigree} \
  --region ${region} \
  --nbins 20 \
  --output ${out_file}


