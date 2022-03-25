#!/usr/bin/env bash
#
#$ -N array_filter
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/array_filter.log
#$ -e logs/array_filter.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 1-9
#$ -tc 1
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"
readonly bash_script="scripts/permute/_array_filter.sh"

readonly chr="${SGE_TASK_ID}"
readonly in_dir="data/permute/counts"
readonly out_dir="data/permute/genes/chr${chr}"

readonly input_path="${in_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_counts_chr${chr}.mt"
readonly input_type='mt'

readonly maf="maf0to5e-2"
readonly annotation="pLoF_damaging_missense"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_chr${chr}"
readonly out_type="mt"

readonly overview="data/permute/overview/overview.tsv.gz"
readonly nslots=2
readonly queue="short.qc"

readonly n_tasks="$( zcat ${overview} | grep "chr${chr}" | wc -l)"
readonly tasks="1-${n_tasks}"

mkdir -p ${out_dir}

set -x
qsub -N "_chr${chr}_permute" \
    -o "logs/_array_filter.log" \
    -e "logs/_array_filter.errors.log" \
    -q "test.qc" \
    -pe shmem 1 \
    -t ${tasks} \
    -tc 1 \
    "${bash_script}" \
    "${chr}" \
    "${input_path}" \
    "${input_type}" \
    "${out_prefix}" \
    "${out_type}" \
    "${overview}" \
    "${nslots}" \
    "${queue}"
set +x


