#!/usr/bin/env bash
#
#$ -N combine_markers
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/combine_markers.log
#$ -e logs/combine_markers.errors.log
#$ -P lindgren.prjc
#$ -q test.qc
#$ -t 14
#$ -V


set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/04_combine_markers.py"

readonly chr="${SGE_TASK_ID}"
readonly markers_dir="data/conditional/common/spa"
readonly ko_dir="data/knockouts/alt"
readonly out_dir="data/conditional/common/knockout/alt"

readonly input_path="${ko_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense.tsv.gz"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense.tsv.gz"
readonly input_type="vcf"
readonly out_type="vcf"

readonly markers=$(cat ${markers_dir}/*.markers | tr "\n" ",")

mkdir -p ${out_dir}

SECONDS=0
set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --input_path ${input_path} \
   --input_type ${input_type} \
   --markers ${markers} \
   --out_type ${out_type} \
   --out_prefix ${out_prefix} \
   && print_update "Finished filtering imputed genotypes ${out_prefix}" ${SECONDS} \
   || raise_error "Filtering impyted genotypes for for ${out_prefix} failed!"



