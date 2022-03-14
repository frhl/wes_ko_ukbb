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

readonly NUM=${SGE_TASK_ID}
readonly gene="$(zcat ${overview} | grep "chr${chr}" | tail -n +2 | cut -f1 | sed ${NUM}'q;d' )"
readonly permutation="$(zcat ${overview} | grep "chr${chr}" | tail -n +2 | cut -f3 | sed ${NUM}'q;d' )"
readonly n_tasks="$(( ( ${permutation} / ${p_per_job} ) + 1 ))"
readonly tasks="1-${n_tasks}"

if [ ${n_tasks} -le ${max_tasks_allowed} ]; then

  echo "${gene} -> ${permutation} .. submitting ${n_tasks} jobs"
  set -x
  qsub -N "_c${chr}_${gene}" \
      -o "logs/_gene_permute.log" \
      -e "logs/_gene_permute.errors.log" \
      -t ${tasks} \
      -q "short.qc" \
      -pe shmem 3 \
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
  >&2 echo "Too many jobs submitted: ${n_tasks}! Exiting.."
fi





