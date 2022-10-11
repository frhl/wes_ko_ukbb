#!/usr/bin/env bash
#
# @description merge phased vcf with unphased parents for later switch error calculation.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=append_parents
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/append_parents.log
#SBATCH --error=logs/append_parents.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-22

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/hail_utils.sh

readonly hail_script="scripts/phasing/07_append_parents.py"
readonly spark_dir="data/tmp/spark"
readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} )

# fam file for calculating switch errors
readonly pedigree_dir="/well/lindgren/UKBIOBANK/nbaya/resources"
readonly pedigree="${pedigree_dir}/ukb11867_pedigree.fam"
# parental genotypes that were not phased
readonly parents_dir="data/unphased/wes_union_calls"
readonly parents_path="${parents_dir}/ukb_wes_union_calls_200k_chr${chr}_parents.vcf.bgz"
# standard genotypes that were phased
readonly phased_dir="data/phased/wes_union_calls/ligated"
readonly phased_path="${phased_dir}/ukb_eur_wes_union_calls_200k_chr${chr}.vcf.bgz"
# out paths and types
readonly out_dir="data/phased/wes_union_calls/with_parents"
readonly out_prefix="${out_dir}/ukb_eur_wes_union_calls_200k_chr${chr}"
readonly out_vcf="${out_prefix}.vcf.bgz"
readonly out_trio="${out_prefix}.trio"
readonly out_type="vcf"

mkdir -p ${out_dir}

#module load BCFtools/1.12-GCC-10.3.0
#make_tabix "${parents_path}" "tbi"
#bcftools merge ${phased_path} ${parents_path} -O z -o ${out_path}

# combine parents and children in same vcf
# note: we can't combine data using bcftools beacuse some variants
# share the same position (but not same reference alleles)
if [ ! -f "${out_vcf}" ]; then
  module purge
  set_up_hail
  set_up_pythonpath_legacy
  python3 ${hail_script} \
    --parents_path ${parents_path} \
    --phased_path ${phased_path} \
    --out_prefix ${out_prefix} \
    --out_type ${out_type}
fi

# calculate switch errors using trio samples
if [ ! -f "${out_trio}" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix ${out_vcf} "tbi"
  bcftools +trio-switch-rate ${out_vcf} -- -p ${pedigree} > ${out_trio}
  switch_errors_by_site ${out_vcf} ${pedigree}
fi


