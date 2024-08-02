#!/usr/bin/env bash
#
# @description combine whole exome sequences with variants from genotyping array
# @note this method is much faster than doing the combination in Hail. It takes
# around 1h for chromosome 21, wheras hail would take 20-30h.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=wes_union_calls
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/wes_union_calls.log
#SBATCH --error=logs/wes_union_calls.errors.log
#SBATCH --partition=short
#SBATCH --constraint="skl-compat"
#SBATCH --cpus-per-task 1
#SBATCH --array=1-22
#
#$ -N wes_union_calls_bcf
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/wes_union_calls.log
#$ -e logs/wes_union_calls.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qe
#$ -t 5,9
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/prefilter/06_wes_union_calls.py"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_wes_dir="data/unphased/wes/prefilter/200k"
readonly in_wes_file="${in_wes_dir}/ukb_split_wes_200k_chr${chr}_no_parents.vcf.bgz"
readonly in_wes_type="vcf"

readonly in_calls_dir="data/unphased/calls/prefilter_no_maf_cutoff/200k"
readonly in_calls_file="${in_calls_dir}/ukb_split_calls_200k_chr${chr}_no_parents.vcf.bgz"
readonly in_calls_type="vcf"

readonly prefix="ukb_wes_union_calls_chr${chr}"
readonly out_dir="data/unphased/wes_union_calls/prefilter_no_maf_cutoff/200k"

readonly tmp_prefix="${out_dir}/${prefix}_tmp"
readonly tmp_file="${tmp_prefix}.vcf.bgz"
readonly tmp_type="vcf"

readonly sort_file="${out_dir}/${prefix}_sorted.vcf.gz"
readonly out_file="${out_dir}/${prefix}.vcf.gz"

# where to store temporary files during sort
readonly tmp_write_dir="data/tmp/bcf/bcfdx"
# this function is used internally in vcf_concat_sort
readonly bcftools_local="/well/lindgren/flassen/software/samtools/v1.16.1/bcftools-v1.16/bcftools-installed-v1.16/bin/bcftools"

mkdir -p ${out_dir}

if [ ! -f "${out_file}" ]; then
  
  # 1. Sort samples of CALLS file by WES file.
  if [ ! -f "${tmp_file}" ]; then
    set_up_hail
    set_up_pythonpath_legacy
    python3 "${hail_script}" \
       --chrom "${chr}" \
       --input_wes_path "${in_wes_file}" \
       --input_wes_type "${in_wes_type}" \
       --input_calls_path "${in_calls_file}" \
       --input_calls_type "${in_calls_type}" \
       --out_calls_prefix "${tmp_prefix}" \
       --out_calls_type "${tmp_type}"
  else
    >&2 echo "${tmp_file} exists. Skipping."
  fi
  SECONDS=0

  module purge
  module load BCFtools/1.12-GCC-10.3.0
  
  # 2. Tabix the temporary file
  make_tabix "${tmp_file}" "tbi"

  module purge
  module load  HTSlib/1.14-GCC-11.2.0  
  # 3. combine the two tables by concatenating and sorting variants 
  vcf_concat_sort ${tmp_file} ${in_wes_file} ${sort_file} ${tmp_write_dir} \
  && print_update "Finished combining VCFs for chr${chr} using bcftools." "${SECONDS}" \
  || raise_error "$( print_update "Combining VCFs for chr${chr} using bcftools failed." ${SECONDS} )"
  make_tabix ${sort_file}

  # 4. calculate AC/AN (Note this was commented out before)
  bcftools +fill-tags ${sort_file} -Oz -o ${out_file} -- -t AN,AC \
  && print_update "Finished fill tags for chr${chr} using bcftools." "${SECONDS}" \
  || raise_error "$( print_update "Fill tag for chr${chr} using bcftools failed." ${SECONDS} )"
  make_tabix ${out_file}
  
  echo "Success. Writing to ${out_file}.."

  # 5. clean up temporary files
  if [ -f "${out_file}" ]; then
    rm "${sort_file}" 
    rm "${sort_file}.tbi"
    rm "${tmp_file}"
    rm "${tmp_file}.tbi"
  fi

else
  rm -f "${tmp_file}"
  rm -f "${tmp_file}.tbi"
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix ${out_file}
  >&2 echo "Final file (${out_file}) already exists! Skipping.."
fi







