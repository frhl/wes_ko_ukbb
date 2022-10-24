#!/usr/bin/env bash
#
# @description split into chunks of samples that are then pre-phased using whatshap
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=merge_prephased_chunks
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/merge_prephased_chunks.log
#SBATCH --error=logs/merge_prephased_chunks.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=21
#
#$ -N prephase_chunks
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prephase_chunks.log
#$ -e logs/prephase_chunks.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 21
#$ -V

set -o errexit
set -o nounset


source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

# how many samples should there be in each chunk 
readonly samples_per_chunk=100
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

# directories and out paths
readonly out_dir="data/phased/wes_union_calls/prephased/chunks"
readonly out_prefix="${out_dir}/ukb_eur_wes_union_calls_200k_chr${chr}"
readonly merge_list="${out_prefix}_spc${samples_per_chunk}_${queue}.mergelist"
readonly tmp="${out_prefix}.tmp"

# remove duplicates (from debugging the functions)
cat ${merge_list} | sort | uniq  > ${tmp}
# combine VCF
module load BCFtools/1.12-GCC-10.3.0
bcftools merge -l ${tmp} -o "${output_prefix}.vcf"
bgzip "${output_prefix}.vcf"
# clean up temporary files
rm ${tmp}

if [ "${out_type}" == "vcf" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_file}.vcf.gz" "tbi"
fi



