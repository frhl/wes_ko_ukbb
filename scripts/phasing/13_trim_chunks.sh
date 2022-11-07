#!/usr/bin/env bash
#
# @description Prune phased chunk boundaries so that only a few overlapping sites remain.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=trim_chunks
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/trim_chunks.log
#SBATCH --error=logs/trim_chunks.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 6
#SBATCH --array=21
#
#
#$ -N trim_chunks
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/trim_chunks.log
#$ -e logs/trim_chunks.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qe
#$ -t 21
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/phasing/13_trim_chunks.py"
readonly spark_dir="data/tmp/spark"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly main_dir="data/phased/wes_scaffold_calls/200k_from_500k/chunks/shapeit5"
readonly in_dir="${main_dir}/ukb_wes_union_calls_shapeit5_200k_from_500k_chr${chr}-16xshort"
readonly in_prefix_regex="shapeit5_prs100000_pro25000_mprs150000" # need this for regex

readonly out_dir="data/phased/wes_scaffold_calls/200k_from_500k/trimmed"
readonly out_prefix="${out_dir}/ukb_wes_scaffold_calls_200k_from_500k_chr${chr}"
readonly out="${out_prefix}.vcf.bgz"

mkdir -p ${out_dir}

if [ ! -f ${out} ]; then
  set_up_hail
  set_up_pythonpath_legacy
  SECONDS=0
  set -x
  python3 ${hail_script} \
      --in_dir "${in_dir}" \
      --in_ext ".vcf.gz" \
      --in_prefix "${in_prefix_regex}" \
      --new_overlap_size 5000 \
      --out_prefix "${out_prefix}" \
      --out_type "vcf" \
      && print_update "Finished merging phased data for chr${chr}" ${SECONDS} \
      || raise_error "Merging phased data for chr${chr} failed" 
  set +x
else
    print_update "Warning: ${out} already exists! Skipping." | tee /dev/stderr
fi






