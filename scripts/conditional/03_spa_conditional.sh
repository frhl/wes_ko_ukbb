#!/usr/bin/env bash
#
#
#$ -N submit_spa_conditional
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/submit_spa_conditional.log
#$ -e logs/submit_spa_conditional.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -V

set -o errexit
set -o nounset

readonly in_dir="data/conditional/common"
readonly out_dir="data/conditional/common"
readonly pheno_dir="data/phenotypes"
readonly step1_dir="data/saige/output/combined/binary/step1"

readonly pheno_list="${pheno_dir}/UKBB_WES200k_binary_phenotypes_header.txt"
readonly index=${SGE_TASK_ID}
readonly phenotype=$( cut -f${index} ${pheno_list} )

readonly spa_script="scripts/conditional/_spa_conditional.sh"
readonly in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
readonly in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
readonly P_cutoff=0.0001

submit_spa_cond_job() 
{
  set -x
  qsub -N "_spa_cond_${3}" \
    -t 21 \
    -q "short.qc@@short.hge" \
    -pe shmem 1 \
    "${spa_script}" \
    "${in_dir}" \
    "${out_dir}" \
    "${pheno_dir}" \
    "${step1_dir}" \
    "${pheno_list}" \
    "${in_gmat}" \
    "${in_var}" \
    "${1}" \
    "${2}" \
    "${P_cutoff}" \
  set +x
}


annotation="synonymous"
vcf="${in_dir}/211111_intervals_${annotation}_${phenotype}.vcf.bgz"
out_prefix="${out_dir}/211111_spa_conditional_${annotation}_${phenotype}"
submit_spa_cond_job ${vcf} ${out_prefix} ${annotation}





