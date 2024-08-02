#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=switch_pp_subset
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/switch_pp_subset.log
#SBATCH --error=logs/switch_pp_subset.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 4
#SBATCH --array=1
#
#$ -N switch_pp_subset
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/switch_pp_subset.log
#$ -e logs/switch_pp_subset.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -t 1-22
#$ -V

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/hail_utils.sh
source utils/qsub_utils.sh


readonly hail_script_subset="scripts/phasing/phasing/06_switch_pp_subset.py"
readonly hail_script_info="scripts/phasing/phasing/05_parents.py"
readonly spark_dir="data/tmp/spark"
readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

# parameter to filter by
readonly pp_cutoff=0.90

# fam file for calculating switch errors
readonly pedigree_dir="/well/lindgren/UKBIOBANK/nbaya/resources"
readonly pedigree="${pedigree_dir}/ukb11867_pedigree.fam"
# for shapeit5 phasing
readonly children_dir="data/phased/wes_union_calls/200k/shapeit5/ligated"
readonly children_path="${children_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly children_type="vcf"
readonly phased_dir="data/phased/wes_union_calls/200k/shapeit5/parents"
readonly phased_path="${phased_dir}/ukb_wes_union_calls_200k_shapeit5_parents_chr${chr}.vcf.gz"
readonly phased_type="vcf"
readonly out_dir="data/phased/wes_union_calls/200k/shapeit5/switch_pp90"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_shapeit5_chr${chr}_pp90"
# out paths and types
readonly out_vcf="${out_prefix}.vcf.bgz"
readonly out_trio="${out_prefix}.trio"
readonly out_trio_by_site="${out_prefix}.txt"
readonly out_trio_by_site_mac="${out_prefix}.long.mac.new"
readonly out_type="vcf"

mkdir -p ${out_dir}

#if [ ! -f "${out_vcf}" ]; then
#  set_up_hail
#  set_up_pythonpath_legacy
#  python3 ${hail_script_subset} \
#    --phased_path ${phased_path} \
#    --phased_type ${phased_type} \
#    --out_prefix ${out_prefix} \
#    --pp_cutoff ${pp_cutoff}
#fi

# calculate switch errors by site
#if [ ! -f "${out_trio_by_site}" ]; then
#  switch_errors_by_site ${out_vcf} ${pedigree}
#fi

if [ ! -f "${out_trio_by_site_mac}.txt.gz" ]; then
  module purge
  set_up_hail
  set_up_pythonpath_legacy
  python3 ${hail_script_info} \
    --children_path ${children_path} \
    --children_type ${children_type} \
    --trio_path ${out_trio_by_site} \
    --out_prefix ${out_trio_by_site_mac}
fi

# calculate switch errors using trio samples
#if [ ! -f "${out_trio}" ]; then
#  module purge
#  module load BCFtools/1.12-GCC-10.3.0
#  make_tabix ${out_vcf} "tbi"
#  bcftools +trio-switch-rate ${out_vcf} -- -p ${pedigree} > ${out_trio}
#fi




