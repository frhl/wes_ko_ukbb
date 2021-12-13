#!/usr/bin/env bash
#
#$ -N submit_spa_conditional_ko_cts 
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/submit_spa_conditional_ko_cts.log
#$ -e logs/submit_spa_conditional_ko_cts.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qe
#$ -t 37
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly vcf_dir="data/conditional/common/merge_markers"
readonly step1_dir="data/saige/output/combined/cts/step1"
readonly out_dir="data/conditional/common/spa_knockout"
readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"

readonly pheno_list="${pheno_dir}/UKBB_WES200k_cts_phenotypes_header.txt"
readonly spa_script="scripts/_spa_test.sh"

readonly in_prefix="211111"

# select phenotype (1-42)
readonly index=${SGE_TASK_ID}
readonly phenotype=$( cut -f${index} ${pheno_list} )

# setup input phenotypes
readonly in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
readonly in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"

# make directories
mkdir -p ${out_dir}

submit_spa_conditional_job() 
{
  qsub -N "spa_${phenotype}_${category}" \
    -t 9 \
    -q "short.qe" \
    -pe shmem 1 \
    "${spa_script}" \
    "${phenotype}" \
    "${in_vcf}" \
    "${in_vcf}.csi" \
    "${in_gmat}" \
    "${in_var}" \
    "${out_prefix}" \
    "${file_markers}"
}

# submit job with specific consequence
submit_spa_with_csqs()
{
  category=${1?Error: Missing arg1 (consequence)}
  in_vcf="${vcf_dir}/${in_prefix}_${category}_${phenotype}_merged_chrCHR.vcf.bgz"
  out_prefix="${out_dir}/${in_prefix}_${category}_${phenotype}"
  file_markers="${vcf_dir}/${in_prefix}_${category}_${phenotype}_merged_chrCHR.cond_markers"
  print_update "Submitting conditional SPA for ${phenotype} [${category}]."
  submit_spa_conditional_job
}

# submit jobs
#submit_spa_with_csqs "ptv"
submit_spa_with_csqs "ptv_damaging_missense"
#submit_spa_with_csqs "synonymous"
#submit_spa_with_csqs "ptv_ptv_LC"
#submit_spa_with_csqs "ptv_ptv_LC_damaging_missense"

