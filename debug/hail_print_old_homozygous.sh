#!/usr/bin/env bash
#
#$ -N old_homozygous
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/old_homozygous.log
#$ -e logs/old_homozygous.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 20
#$ -q short.qe
#$ -t 21

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories

readonly in_dir_phased="/well/lindgren/UKBIOBANK/nbaya/wes_200k/phase_ukb_wes/data/phased/non_singleton"
#readonly in_dir_phased="data/phased"
readonly in_dir_unphased="data/unphased/post-qc"
readonly vep_dir="data/vep/full"
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/mt"

# hail script
readonly hail_script="debug/hail_print_old_homozygous.py"

# input path
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_phased="${in_dir_phased}/ukb_wes_phased_non_singleton_chr${chr}-24xlong.qc-v4.2.2.vcf.gz"
readonly in_unphased="${in_dir_unphased}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly vep="${vep_dir}/ukb_wes_200k_full_vep_chr${chr}.vcf"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_annotated_chr${chr}"
readonly out="${out_prefix}.mt"

# run hail
set_up_hail
set_up_vep
set_up_pythonpath  
python3 "${hail_script}" \
     --chrom ${chr} \
     --input_phased_path ${in_phased}\
     --input_unphased_path ${in_unphased} \
     --input_phased_type "vcf" \
     --input_unphased_type "mt" \
     --vep_path ${vep} \



print_update "Finished running HAIL for chr${chr}" "${SECONDS}"





