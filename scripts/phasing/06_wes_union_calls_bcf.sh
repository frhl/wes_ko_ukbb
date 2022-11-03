#!/usr/bin/env bash
#
# @description combine whole exome sequences with variants from genotyping array
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=wes_union_calls_bcf
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/wes_union_calls_bcf.log
#SBATCH --error=logs/wes_union_calls_bcf.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=21
#
#$ -N wes_union_calls_bcf
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/wes_union_calls_bcf.log
#$ -e logs/wes_union_calls_bcf.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qe
#$ -t 21
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh


readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/06_wes_union_calls_bcf.py"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_wes_dir="data/unphased/wes/prefilter/new"
readonly in_wes_file="${in_wes_dir}/ukb_wes_prefilter_200k_chr${chr}_no_parents.vcf.bgz"
readonly in_wes_type="vcf"

readonly in_calls_dir="data/unphased/calls/prefilter/200k"
readonly in_calls_file="${in_calls_dir}/ukb_split_calls_200k_chr${chr}_no_parents.vcf.bgz"
readonly in_calls_type="vcf"

readonly prefix="ukb_wes_union_calls_chr${chr}"
readonly out_dir="data/unphased/wes_union_calls/bcftools"
readonly tmp_sorted_samples_prefix="${out_dir}/${prefix}_sorted_by_sample"
readonly tmp_sorted_samples_file="${tmp_sorted_samples_prefix}.vcf.bgz"
readonly tmp_sorted_samples_type="vcf"

readonly tmp_unsorted_variants_prefix="${out_dir}/${prefix}_unsorted_by_varaints"
readonly tmp_unsorted_variants_file="${tmp_unsorted_variants_prefix}.vcf.gz"


readonly out_file="${out_dir}/${prefix}.vcf.gz"


mkdir -p ${out_dir}

# 1. Sort samples in calls file by samples in wes file
if [ ! -f "${tmp_sorted_samples_file}" ]; then
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --chrom "${chr}" \
     --input_wes_path "${in_wes_file}" \
     --input_wes_type "${in_wes_type}" \
     --input_calls_path "${in_calls_file}" \
     --input_calls_type "${in_calls_type}" \
     --out_calls_prefix "${tmp_sorted_samples_prefix}" \
     --out_calls_type "${tmp_sorted_samples_type}"
  make_tabix "${tmp_sorted_samples_file}" "tbi"
else
  >&2 echo "${tmp_sorted_samples_file} exists. Skipping."
fi


module load BCFtools/1.12-GCC-10.3.0

# 2. combine the two vcf files by variamts
if [ ! -f "${tmp_unsorted_variants_file}" ]; then
  >&2 echo "Concatenating.."
  bcftools concat ${tmp_sorted_samples_file} ${in_wes_file} -oZ -o ${tmp_unsorted_variants_file}
fi

# 3. sort the file to avoid error:
# [E::hts_idx_push] Unsorted positions on sequence #1: 46664399 followed by 10414352
# index: failed to create index for "did_this_work.vcf.gz"
if [ ! -f "${out_file}" ]; then
  >&2 echo "Sorting file.."
  bcftools sort ${tmp_unsorted_variants_file} -oZ -o ${out_file}
fi


# 4. index file
if [ ! -f "${out_file}.tbi" ]; then
  >&2 echo "Makign tabix.."
  make_tabix "${out_file}" "tbi"
fi


