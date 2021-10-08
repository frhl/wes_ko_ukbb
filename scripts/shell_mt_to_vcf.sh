#!/usr/bin/env bash
#
#$ -N mt_to_vcf
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/mt_to_vcf.log
#$ -e logs/mt_to_vcf.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 20
#$ -q short.qe
#$ -t 8

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly chr=${SGE_TASK_ID}
# directories
readonly in_dir="/well/lindgren/UKBIOBANK/nbaya/wes_200k/ukb_wes_qc/data/filtered"
readonly out_dir="data/unphased/post-qc"
readonly spark_dir="data/tmp/spark"
readonly in="${in_dir}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly hail_script="utils/subscripts/mt_to_vcf.py"
readonly out_prefix="${out_dir}/ukb_wes_200k_filtered_chr${chr}"

# run hail
set_up_hail
set_up_pythonpath
mkdir -p ${out_dir}
python3 "${hail_script}" \
    --input_path ${in} \
    --out_prefix ${out_prefix}


print_update "Finished running HAIL for chr${chr}" "${SECONDS}"





