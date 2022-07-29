#!/usr/bin/env bash
#
#$ -N _prefilter_phenotypes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_prefilter_phenotypes.log
#$ -e logs/_prefilter_phenotypes.errors.log
#$ -P lindgren.prjc
#$ -q short.qa
#$ -pe shmem 10
#$ -t 20-22
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/conditional/rare/_prefilter_phenotypes.R"

readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly pheno_file=${2?Error: Missing arg2 (pheno_file)}
readonly covar_path=${3?Error: Missing arg3 (covar_path)}
readonly in_vcf=${4?Error: Missing arg4 (covar_path)}
readonly out_prefix=${5?Error: Missing arg5 (out_prefix)}

readonly chr=${SGE_TASK_ID}
readonly in_vcf_chr=$(echo ${in_vcf} | sed -e "s/CHR/${chr}/g")
readonly out_prefix_chr=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")
readonly out_prefix_pheno="${out_prefix_chr}_${phenotype}"

if [ ! -f "${out_prefix_pheno}.txt.gz" ]; then
  set_up_rpy
  set -x
  Rscript "${rscript}" \
   --in_vcf ${in_vcf_chr} \
   --out_prefix ${out_prefix_pheno} \
   --phenotype ${phenotype} \
   --pheno_file ${pheno_file} \
   --covariates ${covar_path}
  set +x
else
  >&2 echo "${out_prefix}.txt.gz already exists. Skipping."
fi



