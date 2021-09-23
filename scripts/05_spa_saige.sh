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
#$ -t 1

module purge
source utils/bash_utils.sh

# directories
readonly vcf_dir="derived/knockouts/all/ptvs_damaging_missense"
readonly var_dir="derived/saige/binary/step1"
readonly out_dir="derived/saige/binary/step2"

# input path
readonly chr=${SGE_TASK_ID}
readonly in_vcf="${vcf_dir}/ukb_wes_200k_maf002_miss005_ptv_chr${chr}_ko.vcf.bgz"
readonly in_csi="${vcf}.csi"
readonly spa_script="utils/_spa_test.sh"

# select phenotype (0-42)
set_up_pythonpath
phenotype=$( python utils/extract_phenos_from_header.py \
    --input ${pheno_file} \
    --index ${SGE_TASK_ID})

# setup input phenotypes
readonly in_gmat="${out_dir}/ukb_wes_200k_${phenotype}.rda"
readonly in_var="${out_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"

# output
readonly out_prefix="${out_dir}/ukb_wes_200k_${phenotype}"

submit_spa_job() {
  qsub -N "_spa_${pheno}" \
    -t 22 \
    -q e \
    -pe shmem 2 \
    ${spa_script} \
    ${phenotype} \
    ${in_vcf} \
    ${in_csi} \
    ${in_gmat} \
    ${in_var} \
    ${out_prefix}
}


