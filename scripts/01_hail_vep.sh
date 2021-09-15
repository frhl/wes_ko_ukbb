#!/usr/bin/env bash
#
# Annotate variants
# Note this link when using pick order: https://www.biostars.org/p/120055/
# Author: Frederik Lassen (2021-06-25)
#
#$ -N hail_vep
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/hail_vep.log
#$ -e logs/hail_vep.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qe
#$ -t 22
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh

# Set variables
readonly chr=${SGE_TASK_ID}
readonly in_dir="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/unphased/unfiltered"
readonly in="${in_dir}/ukb_wes_200k_filtered_chr${chr}.vcf.bgz"
readonly out_dir="data/mt" 
readonly out="${out_dir}/ukb_wes_vep_200_chr${chr}" 
readonly spark_dir="data/tmp/spark"
readonly hail_script='utils/hail_vep_export.sh'

# generate VEP annotations as matrix table
set_up_hail
set_up_vep
set_up_pythonpath
mkdir -p ${out_dir}
python3 "${hail_script}" \
	--input_path ${in} \
    --input_type "vcf" \
    --out_prefix ${out} \
    --out_type "mt"


