#!/usr/bin/env bash
#
#$ -N knockout_l
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/knockout_l.log
#$ -e logs/knockout_l.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 10
#$ -q short.qe
#$ -t 1-22

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir="data/mt"
readonly vep_dir="data/vep/full/"
readonly spark_dir="data/tmp/spark"
readonly out_dir="derived/ptvs_only"

# hail script
readonly hail_script="utils/hail_export.py"

# input path
readonly chr=${SGE_TASK_ID}
readonly in_phased="${in_dir}/ukb_wes_200k_annotated_chr${chr}.mt"
readonly in_unphased="${in_dir}/ukb_wes_200k_annotated_chr${chr}_singletons.mt"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_maf002_miss005_ptv_chr${chr}"
readonly out="${out_prefix}.mt"

# run hail
set_up_hail
set_up_pythonpath
mkdir -p ${out_dir}
python3 "${hail_script}" \
    --chrom ${chr} \
    --input_phased_path ${in_phased}\
    --input_unphased_path ${in_unphased} \
    --input_phased_type "mt" \
    --input_unphased_type "mt" \
    --vep_filter "ptv" \
    --maf_max 0.02 \
    --missing 0.05 \
    --out_prefix ${out_prefix} \
    --export_burden \
    --export_ko_probability \
    --export_ko_rsid \
    --export_fake_vcf

print_update "Finished running HAIL for chr${chr}" "${SECONDS}"




