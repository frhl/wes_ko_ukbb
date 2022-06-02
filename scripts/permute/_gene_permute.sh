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
readonly out_prefix=${3?Error: Missing arg2 (in_vcf)}
readonly out_prefix_success=${4?Error: Missing arg2 (in_vcf)}
readonly seed=${5?Error: Missing arg6 (path prefix for saige output)}
readonly gene=${6?Error: Missing arg6 (path prefix for saige output)}
readonly replicates=${7?Error: Missing arg6 (path prefix for saige output)}
readonly cond_markers=${8?Error: Missing arg6 (path prefix for saige output)}

readonly id=${SGE_TASK_ID}
readonly sge_seed=$(( ${id} * ${seed}))
readonly out_prefix_id="${out_prefix}_${id}"
readonly out_file_success="${out_prefix_success}_${id}.SUCCESS"

if [ -f "${input_path}" ]; then
  if [ ! -f "${out_prefix_id}.vcf.gz" ]; then
    set_up_rpy
    Rscript ${rscript} \
      --chrom "chr${chr}" \
      --input_path ${input_path} \
      --input_path_cond ${cond_markers} \
      --permutations ${replicates} \
      --out_prefix ${out_prefix_id} \
      --vcf_id ${gene} \
      --seed ${sge_seed} \
      && print_update "Finished permuting phase for chr${chr}" ${SECONDS} \
      || raise_error "Permuting phase for chr${chr} failed"
    module purge
    module load BCFtools/1.12-GCC-10.3.0
    bgzip "${out_prefix_id}.vcf"
    rm -f "${out_prefix_id}.vcf"
    make_tabix "${out_prefix_id}.vcf.gz" "csi"
  else
    >&2 echo "Error: ${out_prefix_id}.vcf.bgz already exists. Skipping.."
  fi
else
  >&2 echo "Error: ${input_path} does not exist!"
fi


