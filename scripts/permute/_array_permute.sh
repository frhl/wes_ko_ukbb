#!/usr/bin/env bash
#
#
#$ -N _array_permute
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_array_permute.log
#$ -e logs/_array_permute.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe
#$ -V

set -o errexit
set -o nounset

readonly bash_script="scripts/permute/_gene_permute.sh"

readonly chr=${1?Error: Missing arg1 (phenotype)}
readonly input_path=${2?Error: Missing arg1 (phenotype)}
readonly input_type=${3?Error: Missing arg1 (phenotype)}
readonly out_prefix=${4?Error: Missing arg2 (in_vcf)}
readonly out_type=${5?Error: Missing arg3 (in_csi)}
readonly overview=${6?Error: Missing arg6 (path prefix for saige output)}
readonly seed=${7?Error: Missing arg6 (path prefix for saige output)}
readonly p_per_job=${8?Error: Missing arg6 (path prefix for saige output)}
readonly max_tasks_allowed=${9?Error: Missing arg6 (path prefix for saige output)}
readonly nslots=${10?Error: Missing arg6 (path prefix for saige output)}
readonly queue=${11?Error: Missing arg6 (path prefix for saige output)}

readonly NUM=${SGE_TASK_ID}
readonly gene="$(zcat ${overview} | grep "CH" | grep "chr${chr}" | cut -f1 | sed ${NUM}'q;d' )"
readonly permutation="$(zcat ${overview} | grep "CH" | grep "chr${chr}" | cut -f3 | sed ${NUM}'q;d' )"

#readonly n_tasks="$(( ( ${permutation} / ${p_per_job} ) + 1 ))"
#readonly tasks="1-${n_tasks}"
readonly n_tasks="1"
readonly tasks="1"

if [ ${n_tasks} -le ${max_tasks_allowed} ]; then
  echo "${gene} -> ${permutation} .. submitting ${n_tasks} jobs"
  set -x
  qsub -N "_c${chr}_${gene}" \
      -t ${tasks} \
      -q "${queue}" \
      -pe shmem ${nslots} \
      "${bash_script}" \
      "${chr}" \
      "${input_path}" \
      "${input_type}" \
      "${out_prefix}" \
      "${out_type}" \
      "${seed}" \
      "${gene}" \
      "${p_per_job}" \
      "${n_tasks}"
  set +x
else
  >&2 echo "Error! Tried to submit ${n_tasks} task(s) for ${gene}! Aborting submission.."
fi





