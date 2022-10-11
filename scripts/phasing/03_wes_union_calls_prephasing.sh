#!/usr/bin/env bash
#
# @description phase combined set of UK Biobank whole exome and genotyping array using SHAPEIT4.2
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=wes_union_calls_phasing
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/wes_union_calls_phasing.log
#SBATCH --error=logs/wes_union_calls_phasing.errors.log
#SBATCH --partition=long
#SBATCH --cpus-per-task 24
#SBATCH --array=21

source utils/qsub_utils.sh
source utils/vcf_utils.sh
source utils/bash_utils.sh

readonly in_dir="data/unphased/wes_union_calls"
readonly out_dir="data/phased/wes_union_calls/prephase"

readonly chr="${SLURM_ARRAY_TASK_ID}"
readonly in_file="${in_dir}/ukb_eur_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly out_file="${out_dir}/prephased_ukb_eur_wes_union_calls_200k_chr${chr}.vcf.gz"

readonly grch38="/well/lindgren/flassen/ressources/genome_reference/1kg/GRCh38_full_analysis_set_plus_decoy_hla.fa"
readonly cram_dir="/well/ukbb-wes/cram/oqfe/ukbb-11867"

readonly cram_placeholder="${cram_dir}/SAMPLE_oqfe.cram"


# 
whatshap phase \
  --reference ${grch38} \
  --chromosome chr${chr} \
  -o ${out_file} \
  ${vcf_file} \
  ${cram_files}








