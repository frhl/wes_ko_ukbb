#!/usr/bin/env bash

#SBATCH -A lindgren.prj
#SBATCH -J get_replication_binary_phenos
#SBATCH -D /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus
#SBATCH --output=logs/get_replication_binary_phenos.log
#SBATCH --error=logs/get_replication_binary_phenos.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly r_script="scripts/replication/phenotypes/01_get_replication_binary.R"

readonly pheno_dir="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes"
readonly covar_dir="data/phenotypes"

readonly in_dir="data/phenotypes"
readonly in_file="${in_dir}/bin_matrix_eur.txt"

readonly qced_samples_dir="data/samples"
readonly qced_samples="${qced_samples_dir}/UKB.qced_samples.wes450k.txt"

readonly samples_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/phenotypes"
readonly samples_file="${samples_dir}/dec22_phenotypes_binary_200k.tsv.gz"

readonly out_dir="data/phenotypes"
readonly out_prefix="${out_dir}/qced_bin_matrix_eur.remove.wes200k"
readonly out_prefix_inverted="${out_dir}/qced_bin_matrix_eur.keep.wes200k"

# get list of samples in which we set phenotype to NA
readonly samples_to_redact="${out_dir}/UKB.200k.samples.txt"
zcat ${samples_file} | cut -f1 | tail -n+2 > ${samples_to_redact}

mkdir -p ${out_dir}

set_up_rpy

# remove samples in "samples_to_redact"
Rscript ${r_script} \
  --input_path ${in_file} \
  --qced_samples ${qced_samples} \
  --samples_to_redact ${samples_to_redact} \
  --out_prefix ${out_prefix}

# remove samples that are NOT in "samples_to_redact"
Rscript ${r_script} \
  --input_path ${in_file} \
  --qced_samples ${qced_samples} \
  --samples_to_redact ${samples_to_redact} \
  --out_prefix ${out_prefix_inverted} \
  --invert



