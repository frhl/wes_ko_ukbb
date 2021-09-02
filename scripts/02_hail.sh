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
#$ -t 21

set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh

# directories
readonly in_dir_phased="data/phased"
readonly in_dir_unphased="data/unphased/unfiltered"
readonly vep_dir="data/vep/output"
readonly spark_dir="data/tmp/spark"
readonly out_dir="derived/hail"

# hail script
readonly hail_script="utils/hail_export.py"

# input path
readonly chr=${SGE_TASK_ID}
readonly in_phased="${in_dir_phased}/ukb_wes_200k_phased_chr${chr}.1of1.vcf.gz"
readonly in_unphased="${in_dir_unphased}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly vep="${vep_dir}/ukb_wes_200k_full_vep_chr${chr}.vcf"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_phased_eur_maf002_chr${chr}"
readonly out="${out_prefix}.mt"

set_up_hail
mkdir -p ${out_dir}
python3 "${hail_script}" \
    --chrom ${chr} \
    --input_phased_path ${in_phased}\
    --input_unphased_path ${in_unphased} \
    --input_phased_type "vcf" \
    --input_unphased_type "mt" \
    --get_europeans \
    --maf_max 0.02 \
    --missing 0.05 \
    --out_prefix ${out_prefix} \
    --export_burden


print_update "Finished running HAIL for chr${chr}" "${SECONDS}"





