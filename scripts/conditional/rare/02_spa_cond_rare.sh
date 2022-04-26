#!/usr/bin/env bash
#
#$ -N spa_cond_rare
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/spa_cond_rare.log
#$ -e logs/spa_cond_rare.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 1
#$ -tc 10
#$ -V

set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh

readonly vcf_dir="data/knockouts/alt"
readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"

readonly spa_script="scripts/conditional/rare/_spa_cond_rare.sh"
readonly merge_script="scripts/_spa_merge.sh"
readonly in_prefix="ukb_eur_wes_200k"

readonly cond_dir="data/mt/vep"
readonly cond="${cond_dir}/ukb_eur_wes_200k_csqs_chrCHR.tsv.gz"
readonly cond_cat="(pLoF)|(damaging_missense)" # regex

submit_spa_binary_with_csqs()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_spa_with_csqs "${annotation}" "${phenotype}" "binary"
}

submit_spa_cts_with_csqs()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_spa_with_csqs "${annotation}" "${phenotype}" "cts"
}

submit_spa_with_csqs()
{
  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local trait=${3?Error: Missing arg3 (trait)}
  if [ ! -z ${phenotype} ]; then

    local step1_dir="data/saige/output/${trait}/step1"
    local step2_dir="data/saige/output/${trait}/step2"
    local in_vcf="${vcf_dir}/${in_prefix}_chrCHR_${maf}_${annotation}.vcf.bgz"

    if [ "${use_prs}" -eq "0" ]; then
      local in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
      local in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
      local out_prefix="${step2_dir}/${in_prefix}_chrCHR_${maf}_${phenotype}_${annotation}_rare_cond"
      local out_mrg="${step2_dir}/${in_prefix}_${maf}_${phenotype}_${annotation}_rare_cond.txt.gz"
    else
      local in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}_chrCHR.rda"
      local in_var="${step1_dir}/ukb_wes_200k_${phenotype}_chrCHR.varianceRatio.txt"
      local out_prefix="${step2_dir}/${in_prefix}_chrCHR_${maf}_${phenotype}_${annotation}_rare_cond_locoprs"
      local out_mrg="${step2_dir}/${in_prefix}_${maf}_${phenotype}_${annotation}_rare_cond_locoprs.txt.gz"
    fi

    if [ ! -f "${out_mrg}" ]; then
      local qsub_spa_name="spa_${phenotype}_${annotation}"
      local qsub_merge_name="_mrg_${phenotype}_${annotation}"
      submit_spa_job
      submit_merge_job
    else
      >&2 echo "Phenotype ${phenotype} with annotation ${annotation} already exists! Skipping.." 
    fi
  else
    >&2 echo "No phenotype at index ${SGE_TASK_ID}. Exiting.." 
  fi 
}


submit_spa_job() {
  mkdir -p ${step2_dir}
  set -x
  qsub -N "${qsub_spa_name}" \
    -t ${tasks} \
    -q "${queue}" \
    -tc 11 \
    -pe shmem ${nslots} \
    "${spa_script}" \
    "${phenotype}" \
    "${in_vcf}" \
    "${in_vcf}.csi" \
    "${in_gmat}" \
    "${in_var}" \
    "${min_mac}" \
    "${out_prefix}" \
    "${cond}" \
    "${cond_cat}"
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
    "${out_mrg}" \
    "${remove_by_chr}"
  set +x

}

# parameters
readonly conditioning_markers=""
readonly use_prs="1"
readonly min_mac=4
readonly tasks=21
readonly queue="short.qe"
readonly nslots=1



# Binary traits
maf="maf0to5e-2"

# cts traits
submit_spa_cts_with_csqs "pLoF_damaging_missense"
#submit_spa_binary_with_csqs "pLoF_damaging_missense"

#sleep 10
#submit_spa_cts_with_csqs "pLoF"
#submit_spa_binary_with_csqs "pLoF"

#sleep 10
#submit_spa_cts_with_csqs "damaging_missense"
#submit_spa_binary_with_csqs "damaging_missense"

#sleep 10
#submit_spa_cts_with_csqs "synonymous"
#submit_spa_binary_with_csqs "synonymous"





