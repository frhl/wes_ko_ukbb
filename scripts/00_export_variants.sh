#!/usr/bin/env bash
#
# Author: Frederik Lassen (2021-06-25)
#
#$ -N vep_rows
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/vep_rows.log
#$ -e logs/vep_rows.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 10
#$ -q short.qe
#$ -t 1-21
#$ -V


source utils/qsub_utils.sh
source utils/hail_utils.sh

# Set variables
readonly chr=${SGE_TASK_ID}
readonly in_dir="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/unphased/unfiltered"
readonly in="${in_dir}/ukb_wes_200k_filtered_chr${chr}.vcf.bgz"
readonly out_dir="data/variants" 
readonly out="${out_dir}/ukb_wes_200_chr${chr}" 
readonly spark_dir="data/tmp/spark"
readonly hail_script='utils/hail_export_rows.sh'
readonly vep_dir="data/vep/full/"
readonly vep="${vep_dir}/ukb_wes_200k_full_vep_chr${chr}.vcf"

# generate VEP annotations as matrix table
set_up_hail
set_up_vep
set_up_pythonpath
mkdir -p ${out_dir}
python3 "${hail_script}" \
	--input_path ${in} \
    --input_type "vcf" \
    --vep_path "${vep}" \
    --out_prefix ${out}


