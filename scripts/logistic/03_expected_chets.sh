#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=expected_chets
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/expected_chets.log
#SBATCH --error=logs/expected_chets.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1


set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/bash_utils.sh

readonly rscript="scripts/logistic/03_expected_chets.R"

readonly out_dir="data/fisher/"
readonly out_prefix="${out_dir}/expected_chets_chet_hom_eur_wes_200k_pLoF_damaging_missense"
readonly path_phenotypes="data/phenotypes/dec22_phenotypes_binary_200k.tsv.gz"

# parameters
readonly variant_annotation="pLoF_damaging_missense"
# note that "_" charactters are subbed to " "
readonly ko_definition="Homozygote,Compound_heterozygote"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --path_phenotypes ${path_phenotypes} \
  --ko_definition ${ko_definition} \
  --variant_annotation ${variant_annotation} \
  --out_prefix ${out_prefix} 

