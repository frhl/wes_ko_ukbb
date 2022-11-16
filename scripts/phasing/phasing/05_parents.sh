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
#SBATCH --partition=long
#SBATCH --cpus-per-task 2
#SBATCH --array=1,2
#
#$ -N parents
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/parents.log
#$ -e logs/parents.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q long.qc
#$ -t 11
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
readonly parents_dir="data/unphased/wes_union_calls/prefilter_no_maf_cutoff/200k"
readonly parents_path="${parents_dir}/ukb_wes_union_calls_chr${chr}_parents.vcf.gz"
# standard genotypes that were phased
readonly phased_dir="data/phased/wes_union_calls/200k/shapeit5/ligated"
readonly phased_path="${phased_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.bgz"

#readonly phased_dir="data/phased/wes_scaffold_calls/200k_from_500k/ligated"
#readonly phased_path="${phased_dir}/ukb_wes_scaffold_calls_200k_from_500k_chr${chr}.vcf.bgz"
# out paths and types
readonly out_dir="data/phased/wes_union_calls/200k/shapeit5/parents/long-queue"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_shapeit5_parents_chr${chr}"
readonly out_vcf="${out_prefix}.vcf.gz"
readonly out_trio="${out_prefix}.trio"
readonly out_type="vcf"

mkdir -p ${out_dir}

module load BCFtools/1.12-GCC-10.3.0
make_tabix "${parents_path}" "tbi"
make_tabix "${phased_path}" "tbi"
bcftools merge ${phased_path} ${parents_path} -Oz -o ${out_vcf}

# combine parents and children in same vcf
# note: we can't combine data using bcftools beacuse some variants
# share the same position (but not same reference alleles)
#if [ ! -f "${out_vcf}" ]; then
#  module purge
#  set_up_hail
#  set_up_pythonpath_legacy
#  python3 ${hail_script} \
#    --parents_path ${parents_path} \
#    --phased_path ${phased_path} \
#    --out_prefix ${out_prefix} \
#    --out_type ${out_type}
#fi

# calculate switch errors using trio samples
if [ ! -f "${out_trio}" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix ${out_vcf} "tbi"
  bcftools +trio-switch-rate ${out_vcf} -- -p ${pedigree} > ${out_trio}
  switch_errors_by_site ${out_vcf} ${pedigree}
fi


