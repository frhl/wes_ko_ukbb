#!/usr/bin/env bash
#
#$ -N make_genes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/make_genes.log
#$ -e logs/make_genes.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 21
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"
readonly bash_script="scripts/permute/_make_genes.sh"

readonly chr="${SGE_TASK_ID}"
readonly in_dir="data/permute/counts"
readonly out_dir="data/permute/genes/chr${chr}"

readonly input_path="${in_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_counts_chr${chr}.mt"
readonly input_type='mt'

readonly maf="maf0to5e-2"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_chr${chr}"
readonly out_type="mt"

readonly overview="data/permute/overview/overview_genes.tsv.gz"
readonly n_tasks="$( zcat ${overview} | grep "chr${chr}" | wc -l)"
readonly tasks="1-${n_tasks}"

mkdir -p ${out_dir}

set -x
qsub -N "_chr${chr}_permute" \
    -q "short.qc" \
    -pe shmem 1 \
    -t ${tasks} \
    "${bash_script}" \
    "${chr}" \
    "${input_path}" \
    "${input_type}" \
    "${out_prefix}" \
    "${out_type}" \
    "${overview}"
set +x

