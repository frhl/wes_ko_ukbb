#!/usr/bin/env bash
#
# Run Hail & VEP on phased and unphased data seperately since 
# unnphased data only contains singletons, and phased data only
# contains non-singletons. Generates a seperate matrix-table for
# each.
#
#$ -N hail_vep
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/hail_vep.log
#$ -e logs/hail_vep.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 20
#$ -q short.qf
#$ -t 22

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir_phased="data/phased"
readonly in_dir_unphased="data/unphased/post-qc"
readonly vep_dir="data/vep/full/"
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/mt"

# hail script
readonly hail_script="utils/subscripts/hail_format.py"

# input path
readonly chr=${SGE_TASK_ID}
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
mkdir -p ${out_dir}
python3 "${hail_script}" \
    --chrom ${chr} \
    --input_phased_path ${in_phased}\
    --input_unphased_path ${in_unphased} \
    --input_phased_type "vcf" \
    --input_unphased_type "mt" \
    --vep_path ${vep} \
    --out_prefix ${out_prefix}


print_update "Finished running HAIL for chr${chr}" "${SECONDS}"





