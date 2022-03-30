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

readonly bash_script="scripts/permute/_master_permute.sh"

readonly chr=${1?Error: Missing argX}
readonly input_path=${2?Error: Missing argX}
readonly out_prefix=${3?Error: Missing argX}
readonly pheno_dir=${4?Error: Missing argX}
readonly genes_path=${5?Error: Missing argX}
readonly true_p_path=${6?Error: Missing argX}
readonly n_slots=${7?Error: Missing argX}
readonly n_shuffle=${8?Error: Missing argX}
readonly queue=${9?Error: Missing argX}

readonly index=${SGE_TASK_ID}
readonly gene="$(zcat ${genes_path} | grep "chr${chr}" | cut -f1 | sed ${index}'q;d' )"

# chromosome set before
#readonly input_path_chr=$(echo ${input_path} | sed -e "s/CHR/${chr}/g")
#readonly out_prefix_chr=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")

set -x
qsub -N "_c${chr}_${gene}" \
    -q "${queue}" \
    -pe shmem 1 \
    "${bash_script}" \
    "${chr}" \
    "${input_path_chr}" \
    "${out_prefix_chr}" \
    "${pheno_dir}" \
    "${true_p_path}" \
    "${n_slots}" \
    "${n_shuffle}" \
    "${gene}"
set +x





