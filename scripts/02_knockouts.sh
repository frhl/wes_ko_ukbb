#!/usr/bin/env bash
#
# Combine phased and unphased data to make a probabilistic model
# for human knockouts. This function is versatile and produces
# a many different of summary-level files, including:
#
# 1) A fake VCF with probablistic encoding of KO (for SAIGE+ input)
# 2) A matrix containg the variant consequences in each gene for each sample
# 3) A ko probabiltiy matrix containg individuals, genes and probability of KO
# 4) A ko probability matrix with the above + variants involved in KO.
#
#$ -N knockout
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/knockout.log
#$ -e logs/knockout.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 10
#$ -q short.qe
#$ -t 1-22

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

# directories
readonly in_dir="data/mt"
readonly vep_dir="data/vep/full/"
readonly spark_dir="data/tmp/spark"
readonly out_dir="derived/knockouts/all/210925_synonymous"

# hail script
readonly hail_script="utils/subscripts/hail_knockouts.py"

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
    --vep_filter "synonymous" \
    --maf_max 0.02 \
    --missing 0.05 \
    --out_prefix ${out_prefix} \
    --export_fake_vcf \
    --export_burden \
    --export_ko_probability \
    --export_ko_rsid \

print_update "Finished running HAIL for chr${chr}" "${SECONDS}"

# index with csi
module load BCFtools/1.12-GCC-10.3.0
make_tabix "${out_prefix}_ko.vcf.bgz" "csi"



