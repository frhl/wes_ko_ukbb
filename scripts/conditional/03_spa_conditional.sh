#!/usr/bin/env bash
#
#$ -N spa_conditional
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/spa_conditional.log
#$ -e logs/spa_conditional.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 2-44
#$ -V

set -o errexit
set -o nounset

readonly bash_script="scripts/conditional/_spa_conditional.sh"

readonly interval_dir="data/conditional/common/intervals"
readonly out_dir="data/conditional/common/spa_iter"

readonly pheno_dir="data/phenotypes"
readonly in_prefix="ukb_eur_wes_200k"
readonly maf="0to5e-2"

readonly min_mac=4
readonly max_iter=5
readonly P_cutoff="5e-8" #$( from_sci "5e-8" ) 


submit_binary_analysis()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_cond_spa "${annotation}" "${phenotype}" "binary"
}

submit_cts_analysis()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_cond_spa "${annotation}" "${phenotype}" "cts"
}

submit_cond_spa()
{
  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local trait=${3?Error: Missing arg3 (trait)}

  local step1_dir="data/saige/output/${trait}/step1"
  local step2_dir="data/saige/output/${trait}/step2/minmac${min_mac}"
  local in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
  local in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
  local interval_vcf="${interval_dir}/${in_prefix}_maf${maf}_${phenotype}_${annotation}.vcf.bgz"
  local out_prefix="${out_dir}/${in_prefix}_maf${maf}_${phenotype}_${annotation}_cond"

  echo "interval_vcf: $( ls -l ${interval_vcf})"

  mkdir -p ${out_dir}
  #if [[ "$P_cutoff" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    if [ -f ${interval_vcf} ]; then 
      qsub -N "_cond_${1}" \
        -q "short.qc@@short.hge" \
        -t "${SGE_TASK_ID}" \
        -pe shmem 1 \
        "${bash_script}" \
        "${in_gmat}" \
        "${in_var}" \
        "${interval_vcf}" \
        "${out_prefix}" \
        "${P_cutoff}" \
        "${max_iter}" \
        "${min_mac}"
    else
      >&2 echo "${interval_vcf} (interval) does not exist. Exiting.."
    fi
  #else
  #  >&2 echo "${P_cutoff} is not a valid decimal number"
  #fi

}

submit_cts_analysis "pLoF_damaging_missense"
#submit_binary_analysis "pLoF_damaging_missense"




