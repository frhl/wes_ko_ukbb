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
#readonly in_dir="data/permute/counts"
readonly out_dir="data/permute/test_300/chr${chr}"

readonly input_path="${in_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_chr${chr}_GENE.mt"
#readonly input_path="${in_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_counts_chr${chr}.mt"
readonly input_type='mt'

readonly maf="maf0to5e-2"
readonly annotation="pLoF_damaging_missense"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_permuted_chr${chr}"
readonly out_type="vcf"

readonly overview="data/permute/overview/overview.tsv.gz"
readonly max_allowed_jobs=300
readonly p_per_job=300
readonly seed=134
readonly nslots=1
readonly queue="short.qc"

readonly n_tasks="$( zcat ${overview} | grep "chr${chr}" | wc -l)"
#readonly tasks="1-${n_tasks}"
tasks=1

mkdir -p ${out_dir}

set -x
qsub -N "_chr${chr}_permute" \
    -o "logs/_array_permute.log" \
    -e "logs/_array_permute.errors.log" \
    -q "test.qc" \
    -pe shmem 1 \
    -t ${tasks} \
    "${bash_script}" \
    "${chr}" \
    "${input_path}" \
    "${input_type}" \
    "${out_prefix}" \
    "${out_type}" \
    "${overview}" \
    "${seed}" \
    "${p_per_job}" \
    "${max_allowed_jobs}" \
    "${nslots}" \
    "${queue}"
set +x


