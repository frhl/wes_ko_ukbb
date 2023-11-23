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
#SBATCH --cpus-per-task 3
#SBATCH --array=21
#SBATCH --constraint="hge"
#
#
#$ -N trim_chunks
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/trim_chunks.log
#$ -e logs/trim_chunks.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 3
#$ -q short.qc
#$ -t 20-22
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/phasing/phasing/03_trim_chunks.py"
readonly spark_dir="data/tmp/spark"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )


# For eagle2
#readonly main_dir="data/phased/wes_union_calls/200k/eagle2/chunks"
#readonly in_dir="${main_dir}/ukb_wes_union_calls_eagle2_200k_chr${chr}-20xlong"
#readonly in_prefix_regex="eagle2_prs100000_pro25000_mprs150000" 
#readonly out_dir="data/phased/wes_union_calls/200k/eagle2/trimmed"
#readonly out_prefix="${out_dir}/ukb_wes_union_calls_eagle2_200k_chr${chr}_trim"
#readonly out="${out_prefix}.vcf.bgz"

# For SHAPEIT4
#readonly main_dir="data/phased/wes_union_calls/200k/shapeit4/chunks"
#readonly in_dir="${main_dir}/ukb_wes_union_calls_shapeit4_200k_chr${chr}-20xlong"
#readonly in_prefix_regex="shapeit4_prs100000_pro25000_mprs150000" 
#readonly out_dir="data/phased/wes_union_calls/200k/shapeit4/trimmed"
#readonly out_prefix="${out_dir}/ukb_wes_union_calls_shapeit4_200k_chr${chr}_trim"
#readonly out="${out_prefix}.vcf.bgz"

# For SHAPEIT5
readonly main_dir="data/phased/wes_union_calls/200k/shapeit5/phase_rare_test"
readonly in_dir="${main_dir}/ukb_wes_union_calls_shapeit5_200k_chr${chr}-16xlong"
readonly in_prefix_regex="shapeit5_prs100000_pro25000_mprs150000" # need this for regex
readonly out_dir="data/phased/wes_union_calls/200k/shapeit5/trimmed_test"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_shapeit5_200k_chr${chr}_trim"
readonly out="${out_prefix}.vcf.bgz"

#readonly main_dir="data/phased/wes_scaffold_calls/200k_from_500k/chunks/shapeit5"
#readonly in_dir="${main_dir}/ukb_wes_union_calls_shapeit5_200k_from_500k_chr${chr}-16xshort"
#readonly in_prefix_regex="shapeit5_prs100000_pro25000_mprs150000" # need this for regex
#
#readonly out_dir="data/phased/wes_scaffold_calls/200k_from_500k/trimmed"
#readonly out_prefix="${out_dir}/ukb_wes_scaffold_calls_200k_from_500k_chr${chr}"
#readonly out="${out_prefix}.vcf.bgz"

mkdir -p ${out_dir}

if [ ! -f ${out} ]; then
  echo $PATH
  set_up_hail
  echo $PATH
  set_up_pythonpath_legacy
  echo $PATH
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






