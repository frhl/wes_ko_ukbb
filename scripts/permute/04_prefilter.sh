#!/usr/bin/env bash
#
#$ -N _500_array_permute
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_500_array_permute.log
#$ -e logs/_500_array_permute.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qc
#$ -t 21
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/permute/04_array_permute.py"

readonly in_dir="data/mt/annotated"
readonly out_dir="data/permute/output"
readonly genes="data/permute/overview/overview.tsv.gz"

readonly chr="${SGE_TASK_ID}"
readonly input_file="${in_dir}/ukb_eur_wes_200k_annot_chr${chr}.mt"

readonly maf="maf0to5e-2"
readonly annotation="pLoF_damaging_missense"
readonly out_prefix="${out_dir}/test_ukb_eur_wes_200k_permuted_chr${chr}"

SECONDS=0
set_up_hail
set_up_pythonpath_legacy
python3 ${hail_script} \
  --chrom "chr${chr}" \
  --input_path ${input_file} \
  --input_type "mt" \
  --permutations ${genes} \
  --max_permutations 2000 \
  --out_prefix ${out_prefix} \
  --out_type 'mt' \
  --seed "1995" \
  && print_update "Finished permuting phase for chr${chr}" ${SECONDS} \
  || raise_error "Permuting phase for chr${chr} failed"






