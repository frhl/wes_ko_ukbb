#!/usr/bin/env bash
#
#$ -N annotate_variants_vep_qc
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/annotate_variants_vep_qc.log
#$ -e logs/annotate_variants_vep_qc.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q short.qc@@short.hge
#$ -t 21

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly spark_dir="data/tmp/spark"
readonly in_dir="data/hardcalls"
readonly out_dir="data/variants"
readonly vep_dir="data/vep/full"

readonly samples_dir="data/samples"
readonly variant_dir="data/variants"
readonly sex_dir="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/"

# hail script
readonly hail_script="scripts/QC/05_annotate_variants_vep.py"

# input
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly variants="${variants_dir}/02_prefilter_chr${chr}.keep.variant.ht'"
#readonly samples="${samples_dir}/03_initial_qc.keep.sample_list"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_filtered_chr${chr}"

mkdir -p ${out_dir}
set_up_hail
set_up_vep
set_up_pythonpath  
python3 "${hail_script}" \
     --input_path ${in_file}\
     --input_type "mt" \
     --out_prefix ${out_prefix}\
     --gnomad_path ${gnomad}\
     --vep_path ${vep}


