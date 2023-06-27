#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=extract_p
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/extract_p.log
#SBATCH --error=logs/extract_p.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/permute/05_extract_p.R"

readonly in_dir="data/permute/overview"

# overview of significant gene-phenotypes that have been evaluated
readonly path_trait_genes="${in_dir}/phenotypes_with_2cis_2chets.txt.gz"
# directory to be searched recursively permuted P-values
readonly permute_dir="data/permute/permutations_shuffle"

readonly out_dir="data/permute/combined"
readonly out_prefix="${out_dir}/analysis"

mkdir -p ${out_dir}

set_up_rpy
set -x
Rscript "${rscript}" \
  --path_trait_genes "${path_trait_genes}" \
  --permute_dir "${permute_dir}" \
  --out_prefix "${out_prefix}" 



