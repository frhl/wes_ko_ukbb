#!/usr/bin/env bash
#
# @description Get variant csqs and MAF/MAC count by combined tables
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=export_wes_hets
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/export_wes_hets.log
#SBATCH --error=logs/export_wes_hets.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=20-22
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="scripts/phasing/experimental/04_export_wes_hets.py"

readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} ) 
#readonly in_dir="data/mt/annotated"
readonly in_dir="data/unphased/wes/post-qc"
readonly input_prefix="${in_dir}/ukb_wes_200k_filtered_chr${chr}.mt"

readonly out_dir="data/reads/informative_snvs"
readonly out_prefix="${out_dir}/ukb_wes_200k_filtered_chr${chr}"

readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'
readonly final_variant_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list'

mkdir -p ${out_dir}

set_up_hail
set_up_pythonpath_legacy  
python3 "${hail_script}" \
   --in_file ${input_prefix} \
   --final_sample_list ${final_sample_list} \
   --in_type "mt" \
   --out_prefix ${out_prefix} \
   --out_type "tsv" \
   && print_update "Finished exporting csqs chr${chr}" ${SECONDS} \
   || raise_error "Exporting csqs for chr${chr} failed"




