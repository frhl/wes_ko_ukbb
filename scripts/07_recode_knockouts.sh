#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=recode_knockouts
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/recode_knockouts.log
#SBATCH --error=logs/recode_knockouts.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=22
#SBATCH --requeue
# 
#$ -N recode_knockouts
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/recode_knockouts.log
#$ -e logs/recode_knockouts.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 21
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/07_recode_knockouts.R"

readonly cluster=$( get_current_cluster )
readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly vep_dir="data/mt/vep/worst_csq_by_gene_canonical"
readonly vep_file="${vep_dir}/ukb_eur_wes_union_calls_200k_chr${chr}.pp90.tsv.gz"

readonly in_dir="data/knockouts/alt/pp90/combined"
readonly input_path="${in_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense_all.tsv.gz"

readonly out_dir="data/knockouts/alt/pp90/recoded"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}.pp90.recoded"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --input_path ${input_path} \
  --vep_path ${vep_file} \
  --out_prefix ${out_prefix}



