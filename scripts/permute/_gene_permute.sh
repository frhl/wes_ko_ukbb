#!/usr/bin/env bash
#
#
#$ -N _gene_permute
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_gene_permute.log
#$ -e logs/_gene_permute.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/permute/_gene_permute.py"
readonly rscript="scripts/permute/_gene_permute.R"

readonly chr=${1?Error: Missing arg1 (phenotype)}
readonly input_path=${2?Error: Missing arg1 (phenotype)}
readonly input_type=${3?Error: Missing arg1 (phenotype)}
readonly out_prefix=${4?Error: Missing arg2 (in_vcf)}
readonly out_type=${5?Error: Missing arg3 (in_csi)}
readonly seed=${6?Error: Missing arg6 (path prefix for saige output)}
readonly gene=${7?Error: Missing arg6 (path prefix for saige output)}
readonly replicates=${8?Error: Missing arg6 (path prefix for saige output)}
readonly n_tasks=${9?Error: Missing arg6 (path prefix for saige output)}

readonly id=${SGE_TASK_ID}
readonly sge_seed=$(( ${id} * ${seed}))
readonly out_prefix_gene="${out_prefix}_${gene}_${id}of${n_tasks}"
readonly checkpoint="${out_prefix_gene}_${sge_seed}_checkpoint.mt"
readonly input_path_gene=$(echo ${input_path} | sed -e "s/GENE/${gene}/g")

if [ ! -f "${out_prefix_gene}.vcf.gz" ]; then
  set_up_rpy
  Rscript ${rscript} \
    --chrom "chr${chr}" \
    --input_path ${input_path_gene} \
    --permutations ${replicates} \
    --out_prefix ${out_prefix_gene} \
    --vcf_id ${gene} \
    --seed ${sge_seed} \
    && print_update "Finished permuting phase for chr${chr}" ${SECONDS} \
    || raise_error "Permuting phase for chr${chr} failed"
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  bgzip "${out_prefix_gene}.vcf"
  rm -f "${out_prefix_gene}.vcf"
  make_tabix "${out_prefix_gene}.vcf.gz" "csi"
else
  >&2 echo "Error: ${out_prefix_gene}.vcf.bgz already exists. Skipping.."
fi



