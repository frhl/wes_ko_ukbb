#!/usr/bin/env bash
#
#$ -N array_permute
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/array_permute.log
#$ -e logs/array_permute.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 1
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"

readonly spa_script="scripts/_spa_test.sh"
readonly merge_script="scripts/_spa_merge.sh"
readonly vcf_dir="data/knockouts/alt_filtered"

readonly tasks=1-22
readonly queue="short.qf"
readonly nslots=1

submit_chr_binary()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_chr "${annotation}" "${phenotype}" "binary"
}

submit_chr_cts() 
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_chr "${annotation}" "${phenotype}" "cts"
}

submit_chr()
{

  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local trait=${3?Error: Missing arg3 (trait)}

  local step1_dir="data/saige/output/combined/${trait}/step1"
  local step2_dir="data/saige/output/alt/${trait}/step2/combined"
  local out_prefix="${step2_dir}/${in_prefix}_chrCHR_${maf}_${phenotype}_${annotation}"
  local in_vcf="${vcf_dir}/${in_prefix}_chrCHR_${maf}_${annotation}.vcf.bgz"
  local in_tsv="${vcf_dir}/${in_prefix}_chrCHR_${maf}_${annotation}.tsv.gz"
  local in_spa="${step2_dir}/${in_prefix}_${maf}_${annotation}.txt.gz"
  submit_spa_job
  submit_merge_job
}


submit_chr_job() {
  mkdir -p ${step2_dir}
  set -x
  qsub -N "spa_${phenotype}_${category}" \
    -t ${tasks} \
    -q "${queue}" \
    -pe shmem ${nslots} \
    "${bash_script}" \
    "${in_vcf}" \
    "${in_tsv}" \
    "${in_spa}" \
    "${out_prefix}"
  set +x
}


maf="maf0to5e-2"
submit_spa_binary_with_csqs "pLoF_damaging_missense"

#submit_spa_cts_with_csqs "pLoF_damaging_missense"
#done




