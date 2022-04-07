#!/usr/bin/env bash
#
#$ -N _merge_phenotype
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_merge_phenotype.log
#$ -e logs/_merge_phenotype.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hge
#$ -t 1
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly rscript="scripts/simulation/_merge_phenotype.R"

readonly in_prefix=${1?Error: Missing arg1 (phenotype)}
readonly pheno_file=${2?Error: Missing arg1 (phenotype)}
readonly covariates=${3?Error: Missing arg1 (phenotype)}
readonly out_file=${4?Error: Missing arg3 ()}

SECONDS=0
set_up_rpy
set -x
Rscript "${rscript}" \
   --input_path "${out_prefix}.tsv.gz" \
   --real_phenotype_path "${pheno_file}" \
   --covars_keep "${covariates}" \
   --output_path "${out_file}" \
   && print_update "Finished simulating phenotypes for ${in_prefix}" ${SECONDS} \
   || raise_error "Simulating phenotypes for ${in_prefix} failed"
set +x

