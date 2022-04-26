#!/usr/bin/env bash
#
#$ -N extract_markers
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/extract_markers.log
#$ -e logs/extract_markers.errors.log
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
readonly hail_script="scripts/conditional/common/04_extract_markers.py"

readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'

readonly chr="${SGE_TASK_ID}"
readonly out_dir="data/conditional/common/marker"
readonly out_prefix="${out_dir}/conditional_markers_chr${chr}"
readonly out_type="vcf"

readonly markers_dir="data/conditional/common/spa"
readonly markers=$(cat ${markers_dir}/*.markers | tr "\n" ",")

mkdir -p ${out_dir}

SECONDS=0
set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --chrom ${chr} \
   --markers ${markers} \
   --out_type ${out_type} \
   --out_prefix ${out_prefix} \
   && print_update "Finished filtering imputed genotypes ${out_prefix}" ${SECONDS} \
   || raise_error "Filtering imputed genotypes for for ${out_prefix} failed!"



