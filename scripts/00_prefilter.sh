#!/usr/bin/env bash
#
# @description Annotate main MatrixTables with VEP results
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=prefilter
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/prefilter.log
#SBATCH --error=logs/prefilter.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=1

set -o errexit
set -o nounset

# Source utility scripts
source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

# Define directories and script paths
readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="scripts/00_prefilter.py"

# Retrieve task ID and chromosome number
readonly task_id=$(get_array_task_id)
readonly chr=$(get_chr "${task_id}")

# Define input and output paths
readonly in_dir="data/phased/wes_union_calls/200k/shapeit5/clean_ligated"
readonly input_prefix="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly input_type="vcf"

readonly out_dir="data/mt/prefilter_phased"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}_pp90_maf05.vcf.bgz"
readonly out_type="vcf"

# Get paths for final variants and samples
readonly final_sample_list="data/phenotypes/samples/ukb_wes_ko.imputed.qc.samples"
readonly final_variant_list="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list"

# Path to the file containing variants to exclude
readonly exclude="data/genes/220310_common_plofs_to_exclude.txt"

# Define MAF and PP cutoff thresholds
readonly maf_min=0.00
readonly maf_max=0.05
readonly pp_cutoff=0.90

# Create output directory if it does not exist
mkdir -p "${out_dir}"

# Check if the process has been completed before
if [ ! -f "${out_prefix}.mt/_SUCCESS" ]; then
  
  # set up required paths
  set_up_hail 0.2.97
  set_up_pythonpath_legacy

  # Run the Python script with provided arguments
  python3 "${hail_script}" \
     --input_path "${input_prefix}" \
     --input_type "${input_type}" \
     --out_prefix "${out_prefix}" \
     --out_type "${out_type}" \
     --pp_cutoff "${pp_cutoff}" \
     --maf_min "${maf_min}" \
     --maf_max "${maf_max}" \
     --exclude "${exclude}" \
     --final_sample_list "${final_sample_list}" \
     --final_variant_list "${final_variant_list}"

  # Index the output VCF file
  make_tabix "${out_prefix}.vcf.bgz" "csi"
else
  # Output message if the process has been previously completed
  >&2 echo "${out_prefix}.mt already exists! Skipping"
fi


