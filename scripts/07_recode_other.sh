#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=recode_other
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/recode_other.log
#SBATCH --error=logs/recode_other.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=20-22
# 
#$ -N recode_other
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/recode_other.log
#$ -e logs/recode_other.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1-22
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/07_recode_other.R"

readonly cluster=$( get_current_cluster )
readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly vep_dir="data/mt/prefilter/final_90_loftee"
readonly vep_file="${vep_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.csqs.txt.gz"

readonly csqs="synonymous"
readonly in_dir="data/knockouts/alt/pp90/combined"
readonly input_path="${in_dir}/ukb_eur_wes_200k_chr${chr}_${csqs}_all.tsv.gz"


# csqs is appended to file in R script
readonly out_dir="data/knockouts/alt/pp90/recoded_syn"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}.pp90.recoded"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --input_path ${input_path} \
  --vep_path ${vep_file} \
  --out_prefix ${out_prefix} \
  --csqs ${csqs}



