#!/usr/bin/env bash
#
#
#$ -N _ce
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_array_permute.log
#$ -e logs/_array_permute.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe
#$ -V

set -o errexit
set -o nounset

readonly bash_script="scripts/permute/_gene_spa.sh"

readonly chr=${1?Error: Missing arg1 (phenotype)}
readonly phenotype=${2?Error: Missing arg1 (phenotype)}
readonly in_vcf=${3?Error: Missing arg1 (phenotype)}
readonly in_gmat=${4?Error: Missing arg2 (in_vcf)}
readonly in_var=${5?Error: Missing arg3 (in_csi)}
readonly min_mac=${6?Error: Missing arg3 (in_csi)}
readonly overview=${7?Error: Missing arg6 (path prefix for saige output)}
readonly out_prefix=${8?Error: Missing arg6 (path prefix for saige output)}

readonly p_per_job=400

readonly NUM=${SGE_TASK_ID}
readonly gene="$(zcat ${overview} | grep "chr${chr}" | cut -f1 | sed ${NUM}'q;d' )"
readonly permutation="$(zcat ${overview} | grep "chr${chr}" | cut -f3 | sed ${NUM}'q;d' )"
readonly n_tasks="$(( ( ${permutation} / ${p_per_job} ) + 1 ))"
readonly tasks="1-${n_tasks}"

readonly vcf_gene=$(echo ${in_vcf} | sed -e "s/GENE/${gene}/g")
readonly out_gene=$(echo ${out_prefix} | sed -e "s/GENE/${gene}/g")

echo ${vcf_gene}
echo ${out_gene}
echo ${tasks}

set -x
qsub -N "_c${chr}_${gene}" \
    -o "logs/_gene_spa.log" \
    -e "logs/_gene_spa.errors.log" \
    -t ${tasks} \
    -q "short.qc" \
    -pe shmem 1 \
    "${bash_script}" \
    "${chr}" \
    "${vcf_gene}" \
    "${out_gene}" \
    "${in_gmat}" \
    "${in_var}" \
    "${phenotype}" \
    "${gene}" \
    "${n_tasks}" \
    "${min_mac}" 
set +x




