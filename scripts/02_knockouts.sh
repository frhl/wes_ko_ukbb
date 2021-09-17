#!/usr/bin/env bash
#
#$ -N knockout
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/knockout.log
#$ -e logs/knockout.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qe
#$ -t 21

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir_phased="data/phased"
readonly in_dir_unphased="data/unphased/unfiltered"
readonly vep_dir="data/vep/full/"
readonly spark_dir="data/tmp/spark"
readonly out_dir="derived/tmp"

# hail script
readonly hail_script="utils/hail_export.py"

# input path
readonly chr=${SGE_TASK_ID}
readonly in_phased="${in_dir_phased}/ukb_wes_200k_phased_chr${chr}.1of1.vcf.gz"
readonly in_unphased="${in_dir_unphased}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly vep="${vep_dir}/ukb_wes_200k_full_vep_chr${chr}.vcf"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_phased_all_maf002_no_missing_filter_test_chr${chr}"
readonly out="${out_prefix}.mt"

# run hail
set_up_hail
set_up_pythonpath
mkdir -p ${out_dir}
python3 "${hail_script}" \
    --chrom ${chr} \
    --input_phased_path ${in_phased}\
    --input_unphased_path ${in_unphased} \
    --input_phased_type "vcf" \
    --input_unphased_type "mt" \
    --vep_path ${vep} \
    --vep_filter "damaging_missense" "ptv" \
    --maf_max 0.02 \
    --missing 0.05 \
    --out_prefix ${out_prefix} \
    --export_burden \
    --export_ko_probability \
    --export_fake_vcf


print_update "Finished running HAIL for chr${chr}" "${SECONDS}"





