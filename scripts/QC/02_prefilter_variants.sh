#!/usr/bin/env bash
#
#$ -N prefilter_variants
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prefilter_variants.log
#$ -e logs/prefilter_variants.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -t 1-24

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir="data/hardcalls"
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/variants"

# hail script
readonly hail_script="scripts/QC/02_prefilter_variants.py"

# input path
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_file="${in_dir}/ukb_wes_200k_chr${chr}.mt"

# output path
readonly out_prefix="${out_dir}/02_prefilter_chr${chr}.keep.variant.ht"

mkdir -p ${out_dir}
set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
     --input_path ${in_file}\
     --input_type "mt" \
     --out_prefix ${out_prefix}


