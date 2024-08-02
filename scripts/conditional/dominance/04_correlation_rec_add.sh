#!/usr/bin/env bash
#
# Append pseudo variants with actual variants for downstream conditional analysis.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=correlation_rec_add
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/correlation_rec_add.log
#SBATCH --error=logs/correlation_rec_add.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 3
#SBATCH --array=1-22

# note: chrom1 requires more than 12 epyc cpus?

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly rscript="scripts/conditional/dominance/04_correlation_rec_add.R"

readonly array_idx=$( get_array_task_id )
readonly chr="$( get_chr ${array_idx} )"

readonly in_dir="data/conditional/dominance/combine_encodings"
readonly input_path="${in_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense.vcf.gz"

readonly out_dir="data/conditional/dominance/combine_encodings/"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense.correlation"
readonly out_file="${out_prefix}.txt"

#if [ ! -f "${out_file}" ]; then
  set_up_rpy
  Rscript "${rscript}" \
     --chrom ${chr} \
     --input_path ${input_path} \
     --out_prefix ${out_prefix}
#fi


