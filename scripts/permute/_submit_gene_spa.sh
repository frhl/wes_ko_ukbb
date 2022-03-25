#!/usr/bin/env bash
#
#
#$ -N _submit_gene_spa
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_submit_gene_spa.log
#$ -e logs/_submit_gene_spa.errors.log
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
readonly gene_spa=${8?Error: Missing arg6 (path prefix for saige output)}
readonly out_prefix=${9?Error: Missing arg6 (path prefix for saige output)}

readonly p_per_job=300
readonly NUM=${SGE_TASK_ID}
readonly gene="$(zcat ${overview} | grep "chr${chr}" | cut -f1 | sed ${NUM}'q;d' )"
readonly permutation="$(zcat ${overview} | grep "chr${chr}" | cut -f3 | sed ${NUM}'q;d' )"
readonly n_tasks="$(( ( ${permutation} / ${p_per_job} ) + 1 ))"
readonly tasks="1-${n_tasks}"

readonly vcf_gene=$(echo ${in_vcf} | sed -e "s/GENE/${gene}/g")
readonly out_gene=$(echo ${out_prefix} | sed -e "s/GENE/${gene}/g")

echo "Trying to submit ${gene} with ${vcf_gene} with ${n_tasks} task(s) using ${phenotype}.."

# check whether current gene has been tested in this phenotype
# note, may need to also check for the annotation! 
readonly spa_check="$( zcat ${gene_spa} | grep ${phenotype} | grep ${gene} | wc -l )"
readonly task_limit=2000

if [[ ${n_tasks} -le ${task_limit} ]]; then
  if [[ ${spa_check} -ge 0 ]]; then
    set -x
    qsub -N "_c${chr}_${gene}" \
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
  fi
else
  >&2 echo "[${gene}]: ${tasks} tasks is greater than the limit of ${task_limit}! Submission aborted."
fi



