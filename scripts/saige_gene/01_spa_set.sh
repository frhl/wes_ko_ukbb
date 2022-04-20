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
#$ -tc 1
#$ -V

# 12,18

set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh

#readonly vcf_dir="data/knockouts/alt"
readonly vcf_dir="data/mt/annotated"
readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"

readonly spa_script="scripts/saige_gene/_spa_set.sh"
readonly merge_script="scripts/_spa_merge.sh"
readonly in_prefix="ukb_eur_wes_200k_annot"

readonly group_dir="data/mt/vep"
readonly group="${group_dir}/ukb_eur_wes_200k_csqs_chrCHR.saige"

readonly min_mac=4
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
  local step2_dir="data/saige/output/${trait}/step2_set"
  local in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
  local in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
  local out_prefix="${step2_dir}/${in_prefix}_chrCHR_${phenotype}_${annotation}"
  local out_mrg="${step2_dir}/${in_prefix}_chrCHR_${phenotype}_${annotation}.txt.gz"
  local in_vcf="${vcf_dir}/${in_prefix}_chrCHR.vcf.bgz"
  mkdir -p ${step2_dir}
  if [ ! -f ${out_mrg} ]; then
    local qsub_spa_name="sspa_${phenotype}_${annotation}"
    local qsub_merge_name="_smrg_${phenotype}_${annotation}"  
    submit_spa_set_job
    submit_merge_job
  else
    >&2 echo "${out_mrg} already exists. Skipping.."
  fi
}


submit_spa_set_job() 
{
  set -x
  qsub -N "${qsub_spa_name}" \
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
  local remove_by_chr="Y"
  set -x
  qsub -N "${qsub_merge_name}" \
    -q short.qc@@short.hge \
    -pe shmem 1 \
    -hold_jid "${qsub_spa_name}" \
    "${merge_script}" \
    "${out_prefix}" \
    "${step2_dir}" \
    "${out_mrg}" \
    "${remove_by_chr}"
  set +x

}


readonly maf="maf0to5e-2"

# cts traits
submit_spa_set_cts "pLoF"
#submit_spa_set_binary "pLoF"


