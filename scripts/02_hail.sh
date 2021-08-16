#!/usr/bin/env bash
#
#
#$ -N hail_shell
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/hail_shell.log
#$ -e logs/hail_shell.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qe
#$ -t 22

set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh

SGE_TASK_ID=22

# directories
readonly in_dir="data/phased"
readonly vep_dir="data/vep/output"
readonly spark_dir="data/tmp/spark"
readonly out_dir="derived/hail"

# hail script
readonly hail_script="utils/hail_export.py"

# input path
readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly in="${in_dir}/ukb_wes_200k_phased_chr${chr}.1of1.vcf.gz"
readonly vep="${vep_dir}/ukb_wes_200k_vep_chr${chr}.vcf"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_phased_chr${chr}"
readonly out="${out_prefix}.mt"

# setup hail
set_up_hail

SECONDS=0
mkdir -p ${out_dir}
python3 ${hail_script} \
    --chrom ${chr} \
    --input_path ${in} \
    --input_type "vcf" \
    --vep_path ${vep} \
    --get_europeans \
    --out_prefix ${out_prefix} \
    --vep_variants \
    --ko_matrix \
    --ko_samples 

print_update "Finished running HAIL for chr${chr}" "${SECONDS}"





