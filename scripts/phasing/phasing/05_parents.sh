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
#SBATCH --constraint="skl-compat"
#SBATCH --partition=long
#SBATCH --cpus-per-task 4
#SBATCH --array=20-22
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

readonly hail_script="scripts/phasing/phasing/05_parents.py"

readonly spark_dir="data/tmp/spark"
readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

# fam file for calculating switch errors
readonly pedigree_dir="/well/lindgren/UKBIOBANK/nbaya/resources"
readonly pedigree="${pedigree_dir}/ukb11867_pedigree.fam"
# parental genotypes that were not phased
readonly parents_dir="data/unphased/wes_union_calls/prefilter/200k"
readonly parents_path="${parents_dir}/ukb_wes_union_calls_chr${chr}_parents.vcf.gz"
# fore eagle2 phasing
#readonly phased_dir="data/phased/wes_union_calls/200k/eagle2/ligated"
#readonly phased_path="${phased_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.bgz"
#readonly out_dir="data/phased/wes_union_calls/200k/eagle2/parents"
#readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_eagle2_parents_chr${chr}"
# for shapeit4 phasing
readonly phased_dir="data/phased/wes_union_calls/200k/shapeit4/ligated"
readonly phased_path="${phased_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly out_dir="data/phased/wes_union_calls/200k/shapeit4/parents"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_shapeit4_parents_chr${chr}"
# for shapeit5 phasing
#readonly phased_dir="data/phased/wes_union_calls/200k/shapeit5/ligated"
#readonly phased_path="${phased_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.bgz"
#readonly out_dir="data/phased/wes_union_calls/200k/shapeit5/parents_with_hail_count"
#readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_shapeit5_parents_chr${chr}"
# out paths and types
readonly phased_type="vcf"
readonly out_vcf="${out_prefix}.vcf.gz"
readonly out_trio="${out_prefix}.trio"
readonly out_info="${out_prefix}.info"
readonly out_info_gz="${out_prefix}.info.gz"
readonly out_trio_by_site="${out_prefix}.txt"
readonly out_trio_by_site_mac="${out_prefix}.mac"
readonly out_type="vcf"
mkdir -p ${out_dir}

module purge
module load BCFtools/1.12-GCC-10.3.0

# merging files using BCFtools
if [ ! -f "${out_vcf}" ]; then
  echo "merge parents / children: ${out_vcf}"
  make_tabix "${parents_path}" "tbi"
  make_tabix "${phased_path}" "tbi"
  bcftools merge ${phased_path} ${parents_path} -Oz -o ${out_vcf}
fi

# calculate switch errors by site
if [ ! -f "${out_trio_by_site}" ]; then
  echo "SERs by site: ${out_vcf}"
  switch_errors_by_site ${out_vcf} ${pedigree}
fi

# append MAC/AC count 
if [ ! -f "${out_trio_by_site_mac}.txt.gz" ]; then
  module purge
  set_up_hail
  set_up_pythonpath_legacy
  python3 ${hail_script} \
    --children_path ${phased_path} \
    --children_type ${phased_type} \
    --trio_path ${out_trio_by_site} \
    --out_prefix ${out_trio_by_site_mac}
fi

# calculate switch errors using trio samples (standard non-hacky way)
if [ ! -f "${out_trio}" ]; then
  echo "SERs by trio: ${out_trio}"
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix ${out_vcf} "tbi"
  bcftools +trio-switch-rate ${out_vcf} -- -p ${pedigree} > ${out_trio}
fi






