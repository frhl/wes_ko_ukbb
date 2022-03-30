#!/usr/bin/env bash
#
#$ -N array_permute
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/array_permute.log
#$ -e logs/array_permute.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 21
#$ -tc 1
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"
readonly bash_script="scripts/permute/_array_permute.sh"

readonly chr="${SGE_TASK_ID}"
readonly in_dir="data/permute/genes/chr${chr}"
readonly out_dir="data/permute/permutations/chr${chr}"

readonly input_path="${in_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_chr${chr}_GENE.tsv.gz"

readonly maf="maf0to5e-2"
readonly annotation="pLoF_damaging_missense"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_permuted_chr${chr}_GENE"

readonly genes_path="data/permute/overview/overview_genes.tsv.gz"
readonly true_p_path="data/permute/overview/overview_true_p.tsv.gz"

readonly n_shuffle=100
readonly n_slots=1
readonly queue="short.qe"

readonly n_genes="$( zcat ${genes_path} | grep "chr${chr}" | wc -l)"
readonly sge_tasks="1-${n_genes}"

mkdir -p ${out_dir}

set -x
qsub -N "_chr${chr}_permute" \
    -q "test.qc" \
    -pe shmem 1 \
    -t ${sge_tasks} \
    "${bash_script}" \
    "${chr}" \
    "${input_path}" \
    "${out_prefix}" \
    "${pheno_dir}" \
    "${genes_path}" \
    "${true_p_path}" \
    "${n_slots}" \
    "${n_shuffle}" \
    "${queue}"
set +x


