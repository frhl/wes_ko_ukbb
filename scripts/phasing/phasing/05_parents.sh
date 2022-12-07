#!/usr/bin/env bash
#
# @description merge phased vcf with unphased parents for later switch error calculation.
# @note - takes about ~ 24h with 1 a core
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=parents
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/parents.log
#SBATCH --error=logs/parents.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 2
#SBATCH --array=1
#
#$ -N parents
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/parents.log
#$ -e logs/parents.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -t 1
#$ -V

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/hail_utils.sh
source utils/qsub_utils.sh

#readonly hail_script="scripts/phasing/05_parents.py"
readonly spark_dir="data/tmp/spark"
readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

# fam file for calculating switch errors
readonly pedigree_dir="/well/lindgren/UKBIOBANK/nbaya/resources"
readonly pedigree="${pedigree_dir}/ukb11867_pedigree.fam"
# parental genotypes that were not phased
readonly parents_dir="data/unphased/wes_union_calls/prefilter/200k"
readonly parents_path="${parents_dir}/ukb_wes_union_calls_chr${chr}_parents.vcf.gz"
# for shapeit5 phasing
readonly phased_dir="data/phased/wes_union_calls/200k/shapeit5/ligated"
readonly phased_path="${phased_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly out_dir="data/phased/wes_union_calls/200k/shapeit5/parents_improved"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_shapeit5_parents_chr${chr}"
# for eagle2 and shapeit4 testing
#readonly phased_dir="data/phased/wes_union_calls/200k/shapeit4/ukb_wes_union_calls_shapeit4_200k_chr${chr}-16xshort"
#readonly phased_path="${phased_dir}/shapeit4_prs100000_pro25000_mprs150000.1of1.vcf.gz"
#readonly out_dir="data/phased/wes_union_calls/200k/shapeit4/parents_improved"
#readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_shapeit4_parents_chr${chr}"
# out paths and types
readonly out_tmp_vcf="${out_prefix}_tmp.vcf.gz"
readonly out_vcf="${out_prefix}.vcf.gz"
readonly out_trio="${out_prefix}.trio"
readonly out_info="${out_prefix}.info"
readonly out_info_gz="${out_prefix}.info.gz"
readonly out_trio_by_site="${out_prefix}.txt"
readonly out_type="vcf"

mkdir -p ${out_dir}

module purge
module load BCFtools/1.12-GCC-10.3.0

# annotate original (child) file with fill tags. 
# Note: this step takes ~16+ Hours
if [ ! -f "${out_tmp_vcf}" ]; then
  echo "+fill-tags: ${phased_path}"
  bcftools +fill-tags ${phased_path} -Oz -o ${out_tmp_vcf}
fi

# Now merging vcf using BCFtools instead of hail fro speed
if [ ! -f "${out_vcf}" ]; then
  echo "merge parents / children: ${out_vcf}"
  make_tabix "${parents_path}" "tbi"
  make_tabix "${out_tmp_vcf}" "tbi"
  bcftools merge ${out_tmp_vcf} ${parents_path} -Oz -o ${out_vcf}
fi

# write info AC/AN as seperate file
if [ ! -f "${out_info_gz}" ]; then
  echo "Creating INFO file: ${out_info_gz}"
  zcat ${out_vcf} | grep -v "#" | cut -f3,8 > ${out_info}
  gzip ${out_info}
fi

# calculate switch errors by site
if [ ! -f "${out_trio_by_site}" ]; then
  echo "SERs by site: ${out_vcf}"
  switch_errors_by_site ${out_vcf} ${pedigree}
fi

# calculate switch errors using trio samples
if [ ! -f "${out_trio}" ]; then
  echo "SERs by trio: ${out_trio}"
  make_tabix ${out_vcf} "tbi"
  bcftools +trio-switch-rate ${out_vcf} -- -p ${pedigree} > ${out_trio}
fi






