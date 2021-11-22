#!/usr/bin/env bash
#
#$ -N gene_map
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/gene_map.log
#$ -e logs/gene_map.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc@@short.hge
#$ -t 21

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/genes"

# hail script
readonly hail_script="scripts/summary/01_gene_map.py"

# input paths
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly annotation_table="data/vep/hail/ukb_wes_200k_chr${chr}_vep.ht"

# get final samples / variants from Duncan
readonly final_variant_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list'

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_chr${chr}"

# run hail
mkdir -p ${out_dir}
set_up_hail
set_up_pythonpath_legacy  
python3 "${hail_script}" \
     --input_annotation_path ${annotation_table} \
     --final_variant_list ${final_variant_list}\
     --out_prefix ${out_prefix}

print_update "Finished running HAIL for chr${chr}" "${SECONDS}"





