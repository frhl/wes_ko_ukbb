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
readonly merge_script="scripts/_spa_merge.sh"
#readonly in_prefix="ukb_eur_wes_200k_chrCHR"

readonly min_mac=4

readonly tasks=1-22
readonly queue="short.qf"
readonly nslots=1

submit_spa_binary_with_csqs() {
  
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local step1_dir="data/saige/output/combined/binary/step1"
  local out_dir="data/saige/output/combined/binary/step2"
  local pheno_list="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_spa_with_csqs "${annotation}"
}

submit_spa_cts_with_csqs() {
  
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local step1_dir="data/saige/output/combined/cts/step1"
  local out_dir="data/saige/output/combined/cts/step2"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_spa_with_csqs "${annotation}"
}

submit_spa_with_csqs() {
  
  local category=${1?Error: Missing arg1 (consequence)}
  local in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
  local in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
  local out_prefix="${out_dir}/${in_prefix}_${maf}_${phenotype}_${category}"
  local in_vcf="${vcf_dir}/${in_prefix}_${maf}_${category}.vcf.bgz"
  print_update "Submitting SPA for ${phenotype} [${category}]"
  submit_spa_job
  submit_merge_job
}

submit_spa_job() {
  mkdir -p ${out_dir}
  set -x
  qsub -N "spa_${phenotype}_${category}" \
    -t ${tasks} \
    -q "${queue}" \
    -pe shmem ${nslots} \
    "${spa_script}" \
    "${phenotype}" \
    "${in_vcf}" \
    "${in_vcf}.csi" \
    "${in_gmat}" \
    "${in_var}" \
    "${min_mac}" \
    "${out_prefix}"
  set +x
}


submit_merge_job()
{
  readonly remove_by_chr="Y"
  set -x
  qsub -N "_mrg_${phenotype}" \
    -q short.qc@@short.hge \
    -pe shmem 1 \
    -hold_jid "spa_${phenotype}_${category}" \
    "${merge_script}" \
    "${out_prefix}" \
    "${out_dir}" \
    "${out_prefix}.txt.gz" \
    "${remove_by_chr}"
  set +x

}


# Binary traits

for i in $(seq 2 3); do
  in_dir="data/knockouts/null/iter${i}"
  #in_prefix="${out_dir}/ukb_wes_200k_rand_phase"
  in_prefix="ukb_eur_wes_200k_rand_phase_chrCHR"
  maf="maf0to5e-2"
  submit_spa_binary_with_csqs "pLoF_damaging_missense"
done




