#!/usr/bin/env bash
#
#
#$ -N _gene_filter
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_gene_filter.log
#$ -e logs/_gene_filter.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/permute/_gene_filter.py"

readonly chr=${1?Error: Missing arg1 (phenotype)}
readonly input_path=${2?Error: Missing arg1 (phenotype)}
readonly input_type=${3?Error: Missing arg1 (phenotype)}
readonly out_prefix=${4?Error: Missing arg2 (in_vcf)}
readonly out_type=${5?Error: Missing arg3 (in_csi)}
readonly gene=${6?Error: Missing arg6 (path prefix for saige output)}

readonly out_prefix_gene="${out_prefix}_${gene}"

SECONDS=0
set_up_hail
set_up_pythonpath_legacy
python3 ${hail_script} \
  --input_path ${input_path} \
  --input_type ${input_type} \
  --out_prefix ${out_prefix_gene} \
  --out_type ${out_type} \
  --gene ${gene} \
  && print_update "Finished writing gene for chr${chr}" ${SECONDS} \
  || raise_error "writing gene for chr${chr} failed"

