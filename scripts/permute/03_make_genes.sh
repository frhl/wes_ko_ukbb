#!/usr/bin/env bash
#
#$ -N make_genes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/make_genes.log
#$ -e logs/make_genes.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 1-22
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
readonly out_dir="data/permute/genes/phased_only/chr${chr}"

readonly input_path="${in_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_phased_counts_chr${chr}.mt"
readonly input_type='mt'

readonly maf="maf0to5e-2"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_chr${chr}"
readonly out_type="mt"

readonly min_mac=4
readonly genes="data/permute/overview/min_mac${min_mac}/phased_only/main_genes.tsv.gz"

mkdir -p ${out_dir}

if [ -f "${genes}" ]; then
  readonly n_tasks="$( zcat ${genes} | grep -w "chr${chr}" | wc -l)"
  readonly tasks="1-${n_tasks}"
  qsub -N "_c${chr}_make_genes" \
      -q "short.qc" \
      -pe shmem 1 \
      -t ${tasks} \
      "${bash_script}" \
      "${chr}" \
      "${input_path}" \
      "${input_type}" \
      "${out_prefix}" \
      "${out_type}" \
      "${genes}"
else
  >&2 echo "File ${overview} does not exists! Exiting.."
fi

