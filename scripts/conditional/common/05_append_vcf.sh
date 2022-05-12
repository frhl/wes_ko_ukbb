#!/usr/bin/env bash
#
#$ -N append_vcf
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/append_vcf.log
#$ -e logs/append_vcf.errors.log
#$ -P lindgren.prjc
#$ -q short.qc
#$ -pe shmem 2
#$ -t 14
#$ -V


set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/common/05_append_vcf.py"

readonly chr="${SGE_TASK_ID}"
readonly markers_dir="data/conditional/common/spa"
readonly ko_dir="data/knockouts/alt"
readonly out_dir="data/conditional/common/knockout/alt"

readonly markers_path="${markers_dir}/"
readonly input_path="${ko_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense.vcf.bgz"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense"
readonly markers_type="vcf"
readonly input_type="vcf"
readonly out_type="vcf"

mkdir -p ${out_dir}

SECONDS=0
set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --ko_path ${input_path} \
   --ko_type ${input_type} \
   --mk_path ${markers_path} \
   --mk_type ${markers_type} \
   --out_type ${out_type} \
   --out_prefix ${out_prefix} \
   && print_update "Finished merging knockouts with conditional markers ${out_prefix}" ${SECONDS} \
   || raise_error "Merging conditional markers for ${out_prefix} failed!"


