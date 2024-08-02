#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=count_variants_by_category
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/count_variants_by_category.log
#SBATCH --error=logs/count_variants_by_category.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/variants/02_count_variants_by_category.R"

readonly in_dir="data/vep/counts"
readonly in_file="${in_dir}/UKB.exome_array.varaints.vep95.worst_csq_by_gene_canonical.original.counts.txt.gz"

readonly out_dir="data/vep/counts"
readonly out_prefix="${in_dir}/UKB.exome_array.varaints.vep95.worst_csq_by_gene_canonical.original.counts.summarized"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --in_file ${in_file} \
  --out_prefix "${out_prefix}.protein_coding" \
  --filter_biotype "protein_coding"

#Rscript ${rscript} \
#  --in_file ${in_file} \
 # --out_prefix "${out_prefix}"









