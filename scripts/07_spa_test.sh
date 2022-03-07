#!/usr/bin/env bash
#
#$ -N submit_spa_test
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/submit_spa_test.log
#$ -e logs/submit_spa_test.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 1-50
#$ -V

# 12,18

module purge
source utils/bash_utils.sh

readonly vcf_dir="data/knockouts/alt"
readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"

readonly spa_script="scripts/_spa_test.sh"
readonly in_prefix="ukb_wes_200k_maf00_01"

readonly queue="short.qe"
readonly nslots=3

submit_spa_binary_with_csqs() {
  
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local step1_dir="data/saige/output/combined/binary/step1"
  local out_dir="data/saige/output/combined/binary/step2"
  local pheno_list="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
  local phenotype=$( cut -f${SGE_TASK_ID} ${pheno_list} )
  submit_spa_with_csqs "${annotation}"
}

submit_spa_cts_with_csqs() {
  
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local step1_dir="data/saige/output/combined/cts/step1"
  local out_dir="data/saige/output/combined/cts/step2"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( cut -f${SGE_TASK_ID} ${pheno_list} )
  submit_spa_with_csqs "${annotation}"
}

submit_spa_with_csqs() {
  
  local category=${1?Error: Missing arg1 (consequence)}
  local in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
  local in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
  local out_prefix="${out_dir}/${in_prefix}_${phenotype}_${category}"
  local in_vcf="${vcf_dir}/${in_prefix}_chrCHR_${category}_ko.vcf.bgz"
  print_update "Submitting SPA for ${phenotype} [${category}]"
  submit_spa_job
}

submit_spa_job() {
  mkdir -p ${out_dir}
  set -x
  qsub -N "spa_${phenotype}_${category}" \
    -t 1-22 \
    -q "${queue}" \
    -pe shmem ${nslots} \
    "${spa_script}" \
    "${phenotype}" \
    "${in_vcf}" \
    "${in_vcf}.csi" \
    "${in_gmat}" \
    "${in_var}" \
    "${out_prefix}"
  set +x
}


# Binary traits
submit_spa_binary_with_csqs "pLoF"
#submit_spa_binary_with_csqs "ptv_damaging_missense"
#submit_spa_binary_with_csqs "synonymous"
#submit_spa_binary_with_csqs "ptv_ptv_LC"
#submit_spa_binary_with_csqs "ptv_ptv_LC_damaging_missense"

# cts traits
#submit_spa_cts_with_csqs "ptv"
#submit_spa_cts_with_csqs "ptv_damaging_missense"
#submit_spa_cts_with_csqs "synonymous"
#submit_spa_cts_with_csqs "ptv_ptv_LC"
#submit_spa_cts_with_csqs "ptv_ptv_LC_damaging_missense"


