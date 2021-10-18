#!/usr/bin/env bash
#
#
#$ -N step2_saige
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/step2_saige.log
#$ -e logs/step2_saige.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qe
#$ -t 17-19
#$ -V

module purge
source utils/bash_utils.sh
source utils/hail_utils.sh

# directories
readonly vcf_dir="derived/knockouts/all/210925_ptv_damaging_missense"
readonly step1_dir="derived/saige/binary/step1"
readonly out_dir="derived/saige/binary/step2"
readonly pheno_dir="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/200k"
readonly spark_dir="data/tmp/spark"

# input path
readonly chr=${SGE_TASK_ID}
readonly in_vcf="${vcf_dir}/ukb_wes_200k_maf002_miss005_ptv_chr${chr}_ko.vcf.bgz"
readonly in_csi="${vcf}.csi"
readonly pheno_file="${pheno_dir}/UKBB_WES200k_filtered_binary_phenotypes.tsv"
readonly spa_script="utils/subscripts/_spa_test.sh"
readonly pyscript="utils/subscripts/extract_phenos_from_header.py"

# select phenotype (1-42)
set_up_hail
set_up_pythonpath
phenotype=$( python ${pyscript} \
    --input ${pheno_file} \
    --index ${SGE_TASK_ID})

# setup input phenotypes
readonly in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
readonly in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"

# output
readonly out_prefix="${out_dir}/ukb_wes_200k_${phenotype}"

submit_spa_job() {
  #set -x
  qsub -N "spa_${phenotype}" \
    -t 16 \
    -q "short.qe" \
    -pe shmem 4 \
    ${spa_script} \
    ${phenotype} \
    ${in_vcf} \
    ${in_csi} \
    ${in_gmat} \
    ${in_var} \
    ${out_prefix}
}

submit_spa_job
print_update "Submitted SPA scripts for ${phenotype}"

