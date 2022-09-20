#!/usr/bin/env bash
#
#$ -N export_csqs_by_variant
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/export_csqs_by_variant.log
#$ -e logs/export_csqs_by_variant.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hga
#$ -t 1-21

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="scripts/04_export_csqs.py"
readonly rscript="scripts/04_export_csqs.R"

readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_dir="data/mt/annotated"
readonly input_prefix="${in_dir}/ukb_eur_wes_union_calls_200k_chr${chr}.mt"

readonly out_dir="data/mt/vep/worst_csq_for_variant_canonical"
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
     --by "worst_csq_for_variant_canonical" \
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


