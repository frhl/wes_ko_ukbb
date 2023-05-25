#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=exact_test
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/exact_test.log
#SBATCH --error=logs/exact_test.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1


set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/bash_utils.sh

readonly rscript="scripts/logistic/02_exact_test.R"

readonly out_dir="data/fisher/unrelated"
readonly out_prefix="${out_dir}/exact_test_chet_hom_eur_wes_200k_pLoF_damaging_missense"

readonly path_phenotypes="data/phenotypes/dec22_phenotypes_binary_200k.tsv.gz"
readonly path_unrelated="data/post_hoc/unrelated/ukb_wes_ko_samples.txt.gz"

# parameters
# note that "_" charactters are subbed to " "
readonly ko_definition="Homozygote,Compound_heterozygote"
readonly kos_cutoff=5

mkdir -p ${out_dir}

# create shuffled VCF
set_up_rpy
Rscript ${rscript} \
  --out_prefix ${out_prefix} \
  --path_unrelated ${path_unrelated} \
  --path_phenotypes ${path_phenotypes} \
  --ko_definition ${ko_definition} \
  --kos_cutoff ${kos_cutoff}

