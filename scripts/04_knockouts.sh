#!/usr/bin/env bash
#
# Combine phased and unphased data to make a probabilistic model
# for human knockouts.
#
# 1) A fake VCF with probablistic encoding of KO (for SAIGE+ input)
# 2) A matrix containg the variant consequences in each gene for each sample
# 3) A ko probabiltiy matrix containg individuals, genes and probability of KO
# 4) A ko probability matrix with the above + variants involved in KO.
#
#$ -N knockout
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/submit_knockout.log
#$ -e logs/submit_knockout.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hge
#$ -V

#set -o errexit
#set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

# directories
readonly in_dir="data/mt"
readonly spark_dir="data/tmp/spark"
readonly out_dir="derived/knockouts/211108_test"

# input path
readonly knockout_script="scripts/_knockouts.sh"
readonly in_phased="${in_dir}/ukb_wes_200k_annotated_chrCHR.mt"
readonly in_unphased="${in_dir}/ukb_wes_200k_annotated_chrCHR_singletons.mt"
readonly in_phased_type="mt"
readonly in_unphased_type="mt"

# parameters
readonly af_max=0.02

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_af02_chrCHR"

submit_knockout_job() 
{
  qsub -N "_knockout" \
    -t 21 \
    -q "short.qe" \
    -pe shmem 3 \
    "${knockout_script}" \
    "${in_phased}" \
    "${in_phased_type}" \
    "${in_unphased}" \
    "${in_unphased_type}" \
    "${af_max}" \
    "${1}" \
    "${out_prefix}"
}


# parallize jobs
submit_knockout_job "ptv,damaging_misssense"
#submit_knockout_job "ptv"
#submit_knockout_job "synonymous"

# run hail
#set_up_hail
#set_up_pythonpath_legacy
#mkdir -p ${out_dir}
#python3 "${hail_script}" \
#    --chrom ${chr} \
#    --input_phased_path ${in_phased}\
#    --input_unphased_path ${in_unphased} \
#    --input_phased_type "mt" \
#    --input_unphased_type "mt" \
#    --af_max 0.02 \
#    --missing 0.05 \
#    --use_loftee \
#    --out_prefix ${out_prefix} \
#    --export_saige_vcf \
#    --csqs_category "ptv" "damaging_missense"
#    #--export_ko_probability \
#    #--export_ko_rsid 
#    #--export_saige_vcf
#print_update "Finished running HAIL for chr${chr}" "${SECONDS}"



