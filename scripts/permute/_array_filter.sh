#!/usr/bin/env bash
#
#
#$ -N _chr_permute
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_array_permute.log
#$ -e logs/_array_permute.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe
#$ -V

set -o errexit
set -o nounset

readonly bash_script="scripts/permute/_gene_filter.sh"

readonly chr=${1?Error: Missing arg1 (phenotype)}
readonly input_path=${2?Error: Missing arg1 (phenotype)}
readonly input_type=${3?Error: Missing arg1 (phenotype)}
readonly out_prefix=${4?Error: Missing arg2 (in_vcf)}
readonly out_type=${5?Error: Missing arg3 (in_csi)}
readonly overview=${6?Error: Missing arg6 (path prefix for saige output)}
readonly nslots=${7?Error: Missing arg6 (path prefix for saige output)}
readonly queue=${8?Error: Missing arg6 (path prefix for saige output)}

readonly NUM=${SGE_TASK_ID}
readonly gene="$(zcat ${overview} | grep "chr${chr}" | cut -f1 | sed ${NUM}'q;d' )"

set -x
qsub -N "_c${chr}_${gene}" \
    -o "logs/_gene_filter.log" \
    -e "logs/_gene_filter.errors.log" \
    -t ${NUM} \
    -q "${queue}" \
    -pe shmem ${nslots} \
    "${bash_script}" \
    "${chr}" \
    "${input_path}" \
    "${input_type}" \
    "${out_prefix}" \
    "${out_type}" \
    "${gene}" \
set +x





