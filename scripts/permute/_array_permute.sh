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

readonly chr=${1?Error: Missing arg1}
readonly input_path=${2?Error: Missing arg2}
readonly out_prefix=${3?Error: Missing arg3}
readonly pheno_dir=${4?Error: Missing arg4}
readonly genes_path=${5?Error: Missing arg5}
readonly true_p_path=${6?Error: Missing arg6}
readonly min_mac=${7?Error: Missing arg7}
readonly n_replicates=${8?Error: Missing arg8}
readonly n_start_shuffle=${9?Error: Missing arg9}
readonly n_slots_saige=${10?Error: Missing arg10}
readonly n_slots_permute=${11?Error: Missing arg11}
readonly tick_interval=${12?Error: Missing arg12}
readonly tick_timeout=${13?Error: Missing arg13}
readonly queue_saige=${14?Error: Missing arg14}
readonly queue_permute=${15?Error: Missing arg15}
readonly queue_master=${16?Error: Missing arg16}
readonly assoc_format=${17?Error: Missing arg17}

readonly index=${SGE_TASK_ID}
readonly gene="$(zcat ${genes_path} | grep "chr${chr}" | cut -f1 | sed ${index}'q;d' )"

readonly input_path_gene=$(echo ${input_path} | sed -e "s/GENE/${gene}/g")
readonly out_prefix_gene=$(echo ${out_prefix} | sed -e "s/GENE/${gene}/g")

set -x
qsub -N "_c${chr}_${gene}" \
    -q "${queue_master}" \
    -pe shmem 1 \
    "${bash_script}" \
    "${chr}" \
    "${input_path_gene}" \
    "${out_prefix_gene}" \
    "${pheno_dir}" \
    "${true_p_path}" \
    "${min_mac}" \
    "${n_replicates}" \
    "${n_start_shuffle}" \
    "${n_slots_saige}" \
    "${n_slots_permute}" \
    "${tick_interval}" \
    "${tick_timeout}" \
    "${queue_saige}" \
    "${queue_permute}" \
    "${assoc_format}" \
    "${gene}"
set +x





