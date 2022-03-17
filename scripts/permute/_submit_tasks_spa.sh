#!/usr/bin/env bash
#
#
#$ -N _ce
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_tasks_spa.log
#$ -e logs/_tasks_spa.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe
#$ -V

set -o errexit
set -o nounset

readonly bash_script="scripts/permute/_submit_gene_spa.sh"

readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly in_vcf=${2?Error: Missing arg1 (phenotype)}
readonly in_gmat=${3?Error: Missing arg2 (in_vcf)}
readonly in_var=${4?Error: Missing arg3 (in_csi)}
readonly min_mac=${5?Error: Missing arg3 (in_csi)}
readonly overview=${6?Error: Missing arg6 (path prefix for saige output)}
readonly out_prefix=${7?Error: Missing arg6 (path prefix for saige output)}

readonly chr=${SGE_TASK_ID}
readonly n_tasks="$( zcat ${overview} | grep "chr${chr}" | wc -l)"
readonly tasks="1-${n_tasks}"

readonly vcf_chr=$(echo ${in_vcf} | sed -e "s/CHR/${chr}/g")
readonly out_chr=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")

set -x
qsub -N "_c${chr}_spa" \
    -o "logs/submit_gene_spa.log" \
    -e "logs/submit_gene_spa.errors.log" \
    -t ${tasks} \
    -q "test.qc" \
    -pe shmem 1 \
    "${bash_script}" \
    "${chr}" \
    "${phenotype}" \
    "${vcf_chr}" \
    "${in_gmat}" \
    "${in_var}" \
    "${min_mac}" \
    "${overview}" \
    "${out_chr}"
set +x




