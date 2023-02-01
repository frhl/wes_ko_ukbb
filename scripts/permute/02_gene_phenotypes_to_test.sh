#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=gene_phenotypes_to_test
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/gene_phenotypes_to_test.log
#SBATCH --error=logs/gene_phenotypes_to_test.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/permute/02_gene_phenotypes_to_test.R"

readonly phenotypes="data/phenotypes/dec22_phenotypes_binary_200k.tsv.gz"

readonly in_dir="data/permute/overview"
readonly phenotypes_by_count="${in_dir}/phenotypes_with_5cis_5chets.txt.gz"
readonly genes_by_count="${in_dir}/genes_to_run_5cis_5chets.tsv.gz"

readonly min_chets=2

readonly out_dir="data/permute/overview"
readonly out_prefix="${out_dir}/phenotypes_with_5cis_5chets_${min_chets}chetcases"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --phenotypes_by_count ${phenotypes_by_count} \
  --genes_by_count ${genes_by_count} \
  --phenotypes ${phenotypes} \
  --min_chet_cases ${min_chets} \
  --out_prefix ${out_prefix} \

# create seperate header file containing things to run
zcat "${out_prefix}.tsv.gz" | cut -f1 | tail -n +2 > "${out_prefix}_header.txt"


