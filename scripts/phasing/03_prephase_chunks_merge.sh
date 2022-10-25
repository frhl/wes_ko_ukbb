#!/usr/bin/env bash
#
# @description split into chunks of samples that are then pre-phased using whatshap
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=prephase_chunks_merge
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/prephase_chunks_merge.log
#SBATCH --error=logs/prephase_chunks_merge.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 5
#SBATCH --array=21
#
#$ -N prephase_chunks_merge
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prephase_chunks_merge.log
#$ -e logs/prephase_chunks_merge.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qc
#$ -t 21
#$ -V

set -o errexit
set -o nounset


source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

module load BCFtools/1.12-GCC-10.3.0

# how many samples should there be in each chunk 
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

# parameters relating to file name
readonly queue="short"
readonly samples_per_chunk=100

# directories and out paths
readonly out_dir="data/phased/wes_union_calls/prephased/chunks"
readonly out_prefix="${out_dir}/ukb_eur_wes_union_calls_200k_chr${chr}"
readonly input_list="${out_prefix}_spc${samples_per_chunk}_${queue}.mergelist"
readonly tmp="${input_list}.tmp"

readonly out_vcf="${out_prefix}.vcf"
readonly out_vcf_gz="${out_vcf}.gz"

# remove duplicates (from debugging the functions)
cat ${input_list} | sort | uniq | head -n 1000  > ${tmp}

# combine VCFs fast and make tabix
if [ -f "${tmp}" ]; then
  if [ ! -f "${out_vcf}" ]; then
    bcftools merge -l ${tmp} -o "${out_vcf}"
  fi
else
  raise_error "Merge list '${tmp}' does not exist!"
fi

# bgzip
if [ -f "${out_vcf}" ]; then
  if [ ! -f "${out_vcf_gz}" ]; then
    bgzip "${out_vcf}"
  fi
else
  raise_error "File '${out_vcf}' (.vcf) does not exist!"
fi

# index VCF
if [ -f "${out_vcf_gz}" ]; then
  if [ ! -f "${out_vcf_gz}.tbi" ]; then
    make_tabix "${out_vcf_gz}" "tbi"
  fi
else
  raise_error "File '${out_vcf_gz}' (.vcf.gz) does not exist!"
fi

#mv ${tmp} ${merge_list}

# remove duplicates (from debugging the functions)
#cat ${merge_list} | sort | uniq  > ${tmp}
# combine VCF
#module load BCFtools/1.12-GCC-10.3.0
#bcftools merge -l ${tmp} -o "${output_prefix}.vcf"
#bgzip "${output_prefix}.vcf"
# clean up temporary files

#if [ "${out_type}" == "vcf" ]; then
#  module purge
#  module load BCFtools/1.12-GCC-10.3.0
#  make_tabix "${out_file}.vcf.gz" "tbi"
#fi



