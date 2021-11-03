#!/usr/bin/env bash
#
#
#$ -N prepare_spa_saige
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prepare_spa_saige_binary.log
#$ -e logs/prepare_spa_saige_binary.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qe
#$ -t 10
#$ -V

module purge
source utils/bash_utils.sh
source utils/hail_utils.sh

# directories
readonly vcf_dir="derived/knockouts/211102"
readonly step1_dir="data/saige/output/binary/step1"
readonly out_dir="data/saige/output/binary/step2"
readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"

# input path (note chromosome is substituted in _spa_test.sh)
readonly pheno_list="${pheno_dir}/UKBB_WES200k_binary_phenotypes_header.txt"
readonly spa_script="scripts/_spa_test.sh"
readonly in_prefix="ukb_wes_200k_af50"

# select phenotype (1-42)
readonly index=${SGE_TASK_ID}
readonly phenotype=$( cut -f${index} ${pheno_list} )

# setup input phenotypes
readonly in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
readonly in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"

# output
submit_spa_job() 
{
  qsub -N "spa_${phenotype}_${category}" \
    -t 16 \
    -q "short.qe" \
    -pe shmem 1 \
    "${spa_script}" \
    "${phenotype}" \
    "${in_vcf}" \
    "${in_vcf}.csi" \
    "${in_gmat}" \
    "${in_var}" \
    "${out_prefix}"
}

# PTVs
category="ptv"
out_prefix="${out_dir}/${in_prefix}_${category}_${phenotype}"
in_vcf="${vcf_dir}/${in_prefix}_chrCHR_${category}_ko.vcf.bgz"
submit_spa_job

# PTVs + Damaging Missense
category="ptv_damaging_missense"
out_prefix="${out_dir}/${in_prefix}_${category}_${phenotype}"
in_vcf="${vcf_dir}/${in_prefix}_chrCHR_${category}_ko.vcf.bgz"
submit_spa_job

# Non-coding
category="synonymous"
out_prefix="${out_dir}/${in_prefix}_${category}_${phenotype}"
in_vcf="${vcf_dir}/${in_prefix}_chrCHR_${category}_ko.vcf.bgz"
submit_spa_job

print_update "Submitted SPA scripts for ${phenotype}"

