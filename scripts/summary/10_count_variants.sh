#!/usr/bin/env bash
#
#$ -N count_variants
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/count_variants.log
#$ -e logs/count_variants.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hge
#$ -t 1-22

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir_qc="data/qc_old"
readonly spark_dir="data/tmp/spark"
readonly out_dir="derived/tables/worst_csq_for_variant_canonical"

# hail script
readonly hail_script="scripts/summary/10_count_variants.py"

# input paths
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_qc="${in_dir_qc}/ukb_wes_200k_chr${chr}_variants_unphased.ht"
readonly annotation_table="data/vep/hail/ukb_wes_200k_chr${chr}_vep.ht"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_unphased_chr${chr}"

# run hail
mkdir -p ${out_dir}
set_up_hail
set_up_pythonpath_legacy  
python3 "${hail_script}" \
     --chrom ${chr} \
     --input_qc_path ${in_qc}\
     --input_annotation_path ${annotation_table} \
     --out_prefix ${out_prefix}

print_update "Finished running HAIL for chr${chr}" "${SECONDS}"





