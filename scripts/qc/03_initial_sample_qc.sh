#!/usr/bin/env bash
#
#$ -N initial_sample_qc
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/initial_sample_qc.log
#$ -e logs/initial_sample_qc.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -t 1-24

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir="data/unphased/post-qc"
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/samples"
readonly var_dir="data/variants"

# hail script
readonly hail_script="scripts/QC/03_initial_sample_qc.py"

# input path
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly initial_variants="${var_dir}/02_prefilter_chr${chr}.keep.variant.ht"

# output path
readonly out_prefix="${out_dir}/03_chr${chr}_initial_sample_qc.tsv.bgz"

mkdir -p ${out_dir}
set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
     --input_path ${in_file}\
     --input_type "mt" \
     --out_prefix ${out_prefix}\
     --initial_variant_list ${initial_variants}


