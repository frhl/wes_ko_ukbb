#!/usr/bin/env bash
#
#$ -N annotate_genebass
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/annotate_with_genebass.log
#$ -e logs/annotate_with_genebass.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q short.qc@@short.hge
#$ -t 1-22

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/variants"
readonly in_dir="data/unphased/post-qc"

# hail script
readonly hail_script="scripts/QC/08_annotate_with_genebass.py"

# input path
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly final_variant_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list'
readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_genebass_chr${chr}"

mkdir -p ${out_dir}
set_up_hail
set_up_vep
set_up_pythonpath_legacy  
python3 "${hail_script}" \
     --input_path ${in_file}\
     --input_type "mt" \
     --final_variant_list ${final_variant_list}\
     --final_sample_list ${final_sample_list}\
     --out_prefix ${out_prefix}


