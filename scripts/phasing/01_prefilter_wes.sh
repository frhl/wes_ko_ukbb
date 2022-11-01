#!/usr/bin/env bash
#
# @description filter WES quality-controlled data.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=prefilter_wes
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/prefilter_wes.log
#SBATCH --error=logs/prefilter_wes.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=1-22
#
#$ -N prefilter_wes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prefilter_wes.log
#$ -e logs/prefilter_wes.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -t 1-22
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/01_prefilter_wes.py"
readonly hail_vcf_script="scripts/phasing/_to_mt.py"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_dir="data/unphased/wes/post-qc"
readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly in_type="mt"

readonly out_dir="data/unphased/wes/prefilter/new"
readonly out_prefix="${out_dir}/ukb_wes_prefilter_200k_chr${chr}"
readonly out_type="vcf"

#readonly in_dir="data/unphased/wes/prefilter"
#readonly in_file="${in_dir}/ukb_wes_prefilter_200k_chr${chr}.mt"
#readonly in_type="mt"

#readonly out_dir="data/unphased/wes/prefilter"
#readonly out_prefix="${out_dir}/ukb_wes_prefilter_200k_chr${chr}"
#readonly out_type="vcf"

readonly entry_fields_to_drop="GQ,DP,AD,PL"

# samples overlapping exomes and genotypes
readonly samples_list="data/unphased/overlap/ukb_calls_wes_samples.txt"

mkdir -p ${out_dir}
mkdir -p ${spark_dir}

SECONDS=0
set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --input_path "${in_file}" \
   --input_type "${in_type}" \
   --out_prefix "${out_prefix}" \
   --out_type "${out_type}" \
   --drop_entry_fields "${entry_fields_to_drop}" \
   --extract_samples ${samples_list} \
   --exclude_trio_parents \
   --export_parents \
   --min_mac 1 \
   --missing 0.05

#python3 "${hail_vcf_script}" \
#   --input_path "${in_file}" \
#   --input_type "${in_type}" \
#   --out_prefix "${out_prefix}" \
#   --out_type "${out_type}"

module load BCFtools/1.12-GCC-10.3.0
make_tabix "${out_prefix}.vcf.bgz" "tbi"



