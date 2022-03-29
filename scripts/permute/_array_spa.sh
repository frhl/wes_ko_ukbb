#!/usr/bin/env bash
#
#
#$ -N _array_spa
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_array_spa.log
#$ -e logs/_array_spa.errors.log
#$ -P lindgren.prjc
#$ -q test.qc
#$ -pe shmem 1
#$ -q lindgren.qe
#$ -V

set -o errexit
set -o nounset

readonly bash_script="scripts/permute/_submit_gene_spa.sh"

readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly in_vcf=${2?Error: Missing arg2 (in_vcf)}
readonly in_gmat=${3?Error: Missing arg3 (in_gmat)}
readonly in_var=${4?Error: Missing arg4 (in_var)}
readonly min_mac=${5?Error: Missing arg5 (min_mac)}
readonly overview=${6?Error: Missing arg6 (overview)}
readonly gene_spa=${7?Error: Missing arg7 (gene_spa)}
readonly p_per_job=${8?Error: Missing arg8 (p_per_job)}
readonly out_prefix=${9?Error: Missing arg9 (out_prefix)}
readonly out_dir=${10?Error: Missing arg10 (out_dir)}

readonly chr=${SGE_TASK_ID}
readonly n_tasks="$( zcat ${overview} | grep "CH" | grep "chr${chr}" | wc -l)"
#readonly tasks="1-${n_tasks}"
readonly tasks=15

readonly vcf_chr=$(echo ${in_vcf} | sed -e "s/CHR/${chr}/g")
readonly out_chr=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")
readonly out_dir_chr=$(echo ${out_dir} | sed -e "s/CHR/${chr}/g")

readonly task_limit=1000 # genes

if [[ ${tasks} -le ${task_limit} ]]; then
  mkdir -p ${out_dir_chr}
  set -x
  qsub -N "_c${chr}_spa" \
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
      "${gene_spa}" \
      "${p_per_job}" \
      "${out_chr}"
  set +x
else
  >&2 echo "${tasks} tasks is greater than the limit of ${task_limit}! Submission aborted."
fi


