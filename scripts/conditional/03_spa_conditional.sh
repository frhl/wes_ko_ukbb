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
#$ -t 3
#$ -V

set -o errexit
set -o nounset

readonly in_dir="data/conditional/common/filter_genotypes"
readonly pheno_dir="data/phenotypes"
readonly step1_dir="data/saige/output/combined/binary/step1"
readonly out_dir="data/conditional/common/spa_conditional"

readonly pheno_list="${pheno_dir}/UKBB_WES200k_binary_phenotypes_header.txt"
readonly index=${SGE_TASK_ID}
readonly phenotype=$( cut -f${index} ${pheno_list} )

readonly in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
readonly in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"

readonly max_iter=10
readonly P_cutoff=0.01
readonly spa_script="scripts/conditional/_spa_conditional.sh"

mkdir -p ${out_dir}

submit_spa_cond_job() 
{
  set -x
  qsub -N "_spa_cond_${3}" \
    -q "short.qc@@short.hge" \
    -t "${SGE_TASK_ID}" \
    -pe shmem 1 \
    "${spa_script}" \
    "${in_gmat}" \
    "${in_var}" \
    "${1}" \
    "${2}" \
    "${P_cutoff}" \
    "${max_iter}"
  set +x
}


# submit scripts
#annotation="ptv"
#vcf="${in_dir}/211111_intervals_${annotation}_${phenotype}.vcf.bgz"
#out_prefix="${out_dir}/211111_spa_conditional_${annotation}_${phenotype}"
#submit_spa_cond_job ${vcf} ${out_prefix} ${annotation}

#annotation="ptv_damaging_missense"
#vcf="${in_dir}/211111_intervals_${annotation}_${phenotype}.vcf.bgz"
#out_prefix="${out_dir}/211111_spa_conditional_${annotation}_${phenotype}"
#submit_spa_cond_job ${vcf} ${out_prefix} ${annotation}

annotation="synonymous"
vcf="${in_dir}/211111_intervals_${annotation}_${phenotype}.vcf.bgz"
out_prefix="${out_dir}/211111_TEST_P_spa_conditional_${annotation}_${phenotype}"
submit_spa_cond_job ${vcf} ${out_prefix} ${annotation}





