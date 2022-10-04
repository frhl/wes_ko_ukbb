#!/usr/bin/env bash
#
#
#$ -N _make_genes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_make_genes.log
#$ -e logs/_make_genes.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/permute/_make_genes.py"

readonly chr=${1?Error: Missing arg1 (chr)}
readonly input_path=${2?Error: Missing arg2 (input_path)}
readonly input_type=${3?Error: Missing arg3 (input_type)}
readonly out_prefix=${4?Error: Missing arg4 (out_prefix)}
readonly out_type=${5?Error: Missing arg5 (out_type)}
readonly genes=${6?Error: Missing arg6 (genes)}

readonly NUM=${SLURM_ARRAY_TASK_ID}
readonly gene="$(zcat ${genes} | grep -w "chr${chr}" | grep "ENSG" | cut -f1 | sed ${NUM}'q;d' )"
readonly out_prefix_gene="${out_prefix}_${gene}"
readonly out="${out_prefix}.tsv.gz"

if [ ! -f "${out}" ]; then
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy
  python3 ${hail_script} \
    --chrom ${chr} \
    --input_path ${input_path} \
    --input_type ${input_type} \
    --out_prefix ${out_prefix_gene} \
    --out_type ${out_type} \
    --gene ${gene} \
    && print_update "Finished writing gene for chr${chr}" ${SECONDS} \
    || raise_error "writing gene for chr${chr} failed"
else 
 >&2 echo "${out} already exists. Skipping.."
fi


