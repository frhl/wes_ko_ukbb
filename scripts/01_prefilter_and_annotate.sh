#!/usr/bin/env bash
#
# @description Annotate main MatrixTables with VEP results
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=prefilter_and_annotate
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/prefilter_and_annotate.log
#SBATCH --error=logs/prefilter_and_annotate.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=1-21

# used to be called 01_annotate.sh 
# Note: long.qc@@long.hga with 4 slots required to run full pipeline

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="scripts/01_prefilter_and_annotate.py"

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_dir="data/phased/wes_union_calls/200k/shapeit5/ligated"
readonly input_prefix="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly input_type="vcf"

readonly out_dir="data/mt/annotated/old"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}"

readonly out_type="vcf"
readonly out="${out_prefix}.vcf.gz"

readonly annotation_table="data/vep/hail/ukb_wes_200k_chr${chr}_vep.ht"
readonly final_sample_list='data/phenotypes/samples/ukb_wes_ko.imputed.qc.samples'
readonly final_variant_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list'

mkdir -p ${out_dir}

if [ ! -f ${out} ]; then
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy  
  python3 "${hail_script}" \
     --in_file ${input_prefix}\
     --in_type ${input_type} \
     --input_annotation_path ${annotation_table}\
     --final_sample_list ${final_sample_list} \
     --final_variant_list ${final_variant_list}\
     --out_prefix ${out_prefix} \
     --out_type "${out_type}" \
     --annotate_snp_id
fi





