#!/usr/bin/env bash
#
#$ -N count_variants_by_samples
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/count_variants_by_samples.log
#$ -e logs/count_variants_by_samples.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc@@short.hge
#$ -t 1-22

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir_phased="/well/lindgren/UKBIOBANK/nbaya/wes_200k/phase_ukb_wes/data/phased/non_singleton"
readonly in_dir_unphased="data/unphased/post-qc"
readonly spark_dir="data/tmp/spark"
readonly out_dir="derived/qc"

# hail script
readonly hail_script="scripts/summary/03_count_variants_by_samples.py"

# input paths
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_qc="${in_dir_qc}/ukb_wes_200k_chr${chr}_variants_unphased.ht"
readonly annotation_table="data/vep/hail/ukb_wes_200k_chr${chr}_vep.ht"
readonly in_phased="${in_dir_phased}/ukb_wes_phased_non_singleton_chr${chr}-24xlong.qc-v4.2.2.vcf.gz"
readonly in_unphased="${in_dir_unphased}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'
readonly final_variant_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list'

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_unphased_chr${chr}"

# run hail
mkdir -p ${out_dir}
set_up_hail
set_up_pythonpath_legacy  
python3 "${hail_script}" \
     --chrom ${chr} \
     --input_phased_path ${in_phased}\
     --input_unphased_path ${in_unphased} \
     --input_phased_type "vcf" \
     --input_unphased_type "mt" \
     --final_sample_list ${final_sample_list} \
     --final_variant_list ${final_variant_list} \
     --input_annotation_path ${annotation_table} \
     --out_prefix ${out_prefix}

print_update "Finished running HAIL for chr${chr}" "${SECONDS}"





