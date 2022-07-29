#!/usr/bin/env bash
#
#$ -N _prefilter_merge
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_prefilter_merge.log
#$ -e logs/_prefilter_merge.errors.log
#$ -P lindgren.prjc
#$ -q short.qf
#$ -pe shmem 1
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/conditional/rare/_prefilter_merge.R"

readonly in_prefix=${1?Error: Missing arg1 (in_prefix)}
readonly out_prefix=${2?Error: Missing arg2 (out_prefix)}

if [ -f "${in_prefix}.txt.gz" ]; then
  set_up_rpy
  set -x
  Rscript "${rscript}" \
   --in_vcf ${in_vcf} \
   --out_prefix ${out_prefix_pheno} \
   --phenotype ${phenotype} \
   --pheno_file ${pheno_file} \
   --covariates ${covar_path}
  set +x
else
  >&2 echo "${out_prefix}.txt.gz already exists. Skipping."
fi



