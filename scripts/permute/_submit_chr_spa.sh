#!/usr/bin/env bash
##
#$ -N _submit_chr_spa
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_submit_chr_spa.log
#$ -e logs/_submit_chr_spa.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe
#$ -V

set -o errexit
set -o nounset

readonly bash_script="scripts/permute/_submit_tasks_spa.sh"

readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly in_vcf=${2?Error: Missing arg2 (vcf)}
readonly in_gmat=${3?Error: Missing arg2 (gmat)}
readonly in_var=${4?Error: Missing arg3 (variance ratios)}
readonly min_mac=${5?Error: Missing arg3 (min_mac)}
readonly overview=${6?Error: Missing arg6 (overview)}
readonly out_prefix=${7?Error: Missing arg6 (out_prefix)}

readonly tasks=21

set -x
qsub -N "_submit_tasks_spa" \
    -o "logs/submit_tasks_spa.log" \
    -e "logs/submit_tasks_spa.errors.log" \
    -t ${tasks} \
    -q "test.qc" \
    -pe shmem 1 \
    "${bash_script}" \
    "${phenotype}" \
    "${in_vcf}" \
    "${in_gmat}" \
    "${in_var}" \
    "${min_mac}" \
    "${overview}" \
    "${out_prefix}" 
set +x




