#!/usr/bin/env bash
#
#$ -N phased_prefilter
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prefilter_variants_phased.log
#$ -e logs/prefilter_variants_unphased.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qc
#$ -t 21

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir="data/hardcalls/phased"
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/variants/phased"

# hail script
readonly hail_script="scripts/QC/02_prefilter_variants.py"

# input path
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_file="${in_dir}/ukb_wes_200k_chr${chr}.mt"

# output path
readonly out_prefix="${out_dir}/02_prefilter_chr${chr}.keep.variant.ht"

mkdir -p ${out_dir}
set_up_hail
python3 "${hail_script}" \
     --input_path ${in_file}\
     --input_type "mt" \
     --out_prefix ${out_prefix}


