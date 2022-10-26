#!/usr/bin/env bash
#
# @description generate files of genotyped calls
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=isec_samples
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/isec_samples.log
#SBATCH --error=logs/isec_samples.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/00_isec_samples.py"
readonly chr="21"

# output samples
readonly out_dir="data/unphased/"
readonly out_prefix="${out_dir}/ukb_eur_calls_wes_samples"
# get samples in whole exome seq data
readonly in_wes_dir="data/unphased/wes/post-qc"
readonly in_wes_file="${in_wes_dir}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly in_wes_type="mt"
# get european samples
readonly samples_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --chrom "${chr}" \
   --in_path "${in_wes_file}" \
   --in_type "${in_wes_type}" \
   --out_prefix "${out_prefix}" \
   --extract_samples "${samples_list}" \
   --ancestry "eur" \
   --dataset "calls"



