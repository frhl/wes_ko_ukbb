#!/usr/bin/env bash
#
# @description filter WES quality-controlled data.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=prefilter_wes
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/prefilter_wes.log
#SBATCH --error=logs/prefilter_wes.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 2
#SBATCH --array=1-22
#
#$ -N prefilter_wes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o data/help/logs/prefilter_wes.log
#$ -e data/help/logs/prefilter_wes.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -t 21
#$ -V

module purge
source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/help/tmp/spark"
readonly hail_script="scripts/phasing/prefilter/prefilter.py"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_dir="data/unphased/wes/post-qc"
readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly in_type="mt"

readonly out_dir="data/unphased/wes/prefilter/200k/test"
readonly out_prefix="${out_dir}/ukb_prefilter_wes_200k_chr${chr}"
readonly out_prefix_no_parents="${out_prefix}_no_parents"
readonly out_prefix_parents="${out_prefix}_parents"
readonly out_type="mt"

readonly entry_fields_to_drop="GQ,DP,AD,PL"

# samples that are overlapping WES and CALLS
readonly samples_list="data/unphased/overlap/ukb_calls_wes_samples.txt"

mkdir -p ${out_dir}
mkdir -p ${spark_dir}

if [ ! -f "${out_prefix_no_parents}.vcf.bgz" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --input_path "${in_file}" \
     --input_type "${in_type}" \
     --out_prefix "${out_prefix}" \
     --out_type "${out_type}" \
     --drop_fields "${entry_fields_to_drop}" \
     --extract_samples ${samples_list} \
     --missing 0.05 \
     --min_mac 1
else
  >&2 echo "${out_prefix_no_parents} exists. Skipping."
fi




