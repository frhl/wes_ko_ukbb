#!/usr/bin/env bash
#
#
#$ -N prepare_spa_saige
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prepare_spa_saige_cts.log
#$ -e logs/prepare_spa_saige_cts.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qe
#$ -t 1-45
#$ -V

# use 37/38

module purge
source utils/bash_utils.sh
source utils/hail_utils.sh

# directories
readonly vcf_dir="derived/knockouts/211111"
readonly step1_dir="data/saige/output/combined/cts/step1"
readonly out_dir="data/saige/output/combined/cts/step2/211111"
readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"

# input path (note chromosome is substituted in _spa_test.sh)
readonly pheno_list="${pheno_dir}/UKBB_WES200k_cts_phenotypes_header.txt"
readonly spa_script="scripts/_spa_test.sh"
readonly in_prefix="ukb_wes_200k_maf01_50"

# select phenotype (1-42)
readonly index=${SGE_TASK_ID}
readonly phenotype=$( cut -f${index} ${pheno_list} )

# setup input phenotypes
readonly in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
readonly in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"

# make directories
mkdir -p ${out_dir}

submit_spa_job() 
{
  qsub -N "spa_${phenotype}_${category}" \
    -t 1-22 \
    -q "short.qe" \
    -pe shmem 1 \
    "${spa_script}" \
    "${phenotype}" \
    "${in_vcf}" \
    "${in_vcf}.csi" \
    "${in_gmat}" \
    "${in_var}" \
    "${out_prefix}"
}


# submit job with specific consequence
submit_spa_with_csqs()
{
  category=${1?Error: Missing arg1 (consequence)}
  out_prefix="${out_dir}/${in_prefix}_${category}_${phenotype}"
  in_vcf="${vcf_dir}/${in_prefix}_chrCHR_${category}_ko.vcf.bgz"
  print_update "Submitting SPA for ${phenotype} [${category}]"
  submit_spa_job
}

# submit jobs
submit_spa_with_csqs "ptv"
submit_spa_with_csqs "ptv_damaging_missense"
submit_spa_with_csqs "synonymous"
#submit_spa_with_csqs "ptv_lc"
#submit_spa_with_csqs "ptv_lc_damaging_missense"

