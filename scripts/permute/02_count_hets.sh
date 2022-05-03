#!/usr/bin/env bash
#
#$ -N count_hets
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/count_hets.log
#$ -e logs/count_hets.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -t 1-22
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/permute/02_count_hets.py"

readonly in_dir="data/mt/csqs"
readonly out_dir="data/permute/counts"

readonly chr="${SGE_TASK_ID}"
readonly input_path="${in_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense.mt"
readonly input_type='mt'

readonly annotation="pLoF_damaging_missense"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_counts_chr${chr}"
readonly out_type="mt"

mkdir -p ${out_dir}

if [ ! -d "${out_prefix}.mt" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  python3 ${hail_script} \
    --input_path ${input_path} \
    --input_type ${input_type} \
    --out_prefix ${out_prefix} \
    --out_type ${out_type}
else
  >&2 echo "${out_prefix}.mt already exists. Skipping.."
fi




