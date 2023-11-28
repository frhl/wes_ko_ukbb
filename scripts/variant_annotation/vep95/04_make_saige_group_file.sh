#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=make_saige_group_file
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/make_saige_group_file.log
#SBATCH --error=logs/make_saige_group_file.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/variant_annotation/vep95/04_make_saige_group_file.R"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_dir="data/vep/vep95/worst_csqs"
readonly in_path="${in_dir}/UKB.chr${chr}.exome_array.variants_only.vep95.csqs.worst_csq_by_gene_canonical.original.txt.gz"

readonly out_dir="data/vep/vep95/saige_group"
readonly out_saige="${out_dir}/UKB.chr${chr}.exome_array.variants_only.vep95.csqs.worst_csq_by_gene_canonical.original.saige.txt"

readonly annotations_to_merge="pLoF,damaging_missense"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --input_path "${in_path}" \
  --output_path "${out_saige}" \
  --merge_into_single_annotation "${annotations_to_merge}" \
  --delimiter " "


gzip -f ${out_saige}

