#!/usr/bin/env bash
#
#$ -N export_csqs
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/export_csqs.log
#$ -e logs/export_csqs.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hga
#$ -t 21

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="scripts/03_export_csqs.py"

readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_dir="data/mt/annotated"
readonly input_prefix="${in_dir}/ukb_eur_wes_200k_annot_chr${chr}.mt"

readonly out_dir="data/mt/csqs"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}"


mkdir -p ${out_dir}

set_up_hail
set_up_pythonpath_legacy  
set -x
python3 "${hail_script}" \
   --in_file ${input_prefix}\
   --in_type "mt" \
   --out_prefix ${out_prefix}
set +x 






