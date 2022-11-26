#!/usr/bin/env bash
#
# @description Get variant csqs and MAF/MAC count by combined tables
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=export_csqs
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/export_csqs.log
#SBATCH --error=logs/export_csqs.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=20-22
#SBATCH --requeue
# 
#$ -N export_csqs
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/export_csqs.log
#$ -e logs/export_csqs.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1-19
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="scripts/03_export_csqs.py"
readonly rscript="scripts/03_export_csqs.R"

readonly cluster=$( get_current_cluster )
readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_dir="data/mt/annotated"
readonly input_prefix="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.mt"

readonly out_dir="data/mt/vep/worst_csq_by_gene_canonical"
readonly out_prefix="${out_dir}/ukb_eur_wes_union_calls_200k_chr${chr}"
readonly out_saige="${out_prefix}.saige"


mkdir -p ${out_dir}

if [ ! -f "${out_prefix}.tsv.gz" ]; then
  set_up_hail
  set_up_pythonpath_legacy  
  python3 "${hail_script}" \
     --in_file ${input_prefix}\
     --in_type "mt" \
     --out_prefix ${out_prefix} \
     --by "worst_csq_by_gene_canonical" \
     --by_explode \
     --out_type "tsv" \
     && print_update "Finished exporting csqs chr${chr}" ${SECONDS} \
     || raise_error "Exporting csqs for chr${chr} failed"
else
  >&2 echo "${out_prefix} already exists. Skipping.."
fi 


# Generate SAIGE-GENE+ Group file consequence
# annotations (SAIGE version > 0.99.2)
module purge
set_up_rpy
Rscript ${rscript} \
  --input_path "${out_prefix}.tsv.gz" \
  --output_path "${out_saige}" \
  --delimiter " "


