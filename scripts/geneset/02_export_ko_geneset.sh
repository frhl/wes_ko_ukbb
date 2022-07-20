#!/usr/bin/env bash
#
#$ -N export_ko_geneset
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/export_ko_geneset.log
#$ -e logs/export_ko_geneset.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hga
#$ -t 22

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="scripts/04_export_csqs.py"
readonly rscript="scripts/04_export_csqs_test.R"

readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_dir="data/mt/annotated"
readonly input_prefix="${in_dir}/ukb_eur_wes_200k_annot_chr${chr}.mt"

readonly out_dir="data/mt/vep"
#readonly out_prefix="${out_dir}/ukb_eur_wes_200k_csqs_chr${chr}"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_csqs_chrALL_test"
readonly out_saige="${out_prefix}.saige"


mkdir -p ${out_dir}


# Generate SAIGE-GENE+ Group file consequence
# annotations (SAIGE version > 0.99.2)
module purge
set_up_rpy
Rscript ${rscript} \
  --input_path "${out_prefix}.tsv.gz" \
  --output_path "${out_saige}" \
  --delimiter " "


