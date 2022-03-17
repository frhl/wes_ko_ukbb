#!/usr/bin/env bash
#
#$ -N spa_set
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/spa_set.log
#$ -e logs/spa_set.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 1
#$ -V

# 12,18

module purge
source utils/bash_utils.sh

readonly vcf_dir="data/knockouts/alt"
readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"

readonly spa_script="scripts/_spa_set.sh"
readonly merge_script="scripts/_spa_merge.sh"
readonly in_prefix="ukb_eur_wes_200k_chrCHR"

readonly group_dir="data/mt/vep"
readonly group="${group_dir}/ukb_eur_wes_200k_csqs_chrCHR.saige"

readonly min_mac=5

readonly tasks=21
readonly queue="short.qf"
readonly nslots=1

submit_spa_set_binary()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_spa_pair "${annotation}" "${phenotype}" "binary"
}

submit_spa_set_cts()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_spa_pair "${annotation}" "${phenotype}" "cts"
}

submit_spa_pair()
{

  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local trait=${3?Error: Missing arg3 (trait)}

  local step1_dir="data/saige/output/${trait}/step1"
  local step2_dir="data/saige/output/set/${trait}/step2"
  local in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
  local in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
  local out_prefix="${step2_dir}/${in_prefix}_${maf}_${phenotype}_${annotation}"
  local in_vcf="${vcf_dir}/${in_prefix}_${maf}_${annotation}.vcf.bgz"
  submit_spa_set_job
 # submit_merge_job
}


submit_spa_set_job() 
{
  mkdir -p ${step2_dir}
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
    "${group}" \
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
maf="maf0to5e-2"
#submit_spa_binary_with_csqs "pLoF"
#submit_spa_binary_with_csqs "pLoF_damaging_missense"
#submit_spa_binary_with_csqs "synonymous"

# cts traits
submit_spa_set_cts "pLoF_damaging_missense"
#submit_spa_binary_with_csqs "pLoF_damaging_missense"


