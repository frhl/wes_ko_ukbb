#!/usr/bin/env bash
#
#$ -N aggr_knockouts_by_survival_phenos
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/aggr_knockouts_by_survival_phenos.log
#$ -e logs/aggr_knockouts_by_survival_phenos.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly rscript="scripts/survival/02_aggr_knockouts_by_survival_phenos.R"

readonly pheno_dir="data/phenotypes"
readonly phenos="${pheno_dir}/filtered_phenotypes_binary.tsv"
readonly ko_dir="data/knockouts/alt"
readonly ko_file="${ko_dir}/ukb_eur_wes_200k_chrCHR_maf0to5e-2_pLoF_damaging_missense_all.tsv.gz"
readonly out_dir="data/post_hoc/results"
readonly out_prefix="${out_dir}/pLoF_damaging_missense_maf0to5e-2_survival_knockouts"

readonly samvida_phenos="/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_time_to_event_matrix.txt"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --phenotypes ${samvida_phenos} \
  --original_phenotypes ${phenos} \
  --ko_file ${ko_file} \
  --out_prefix ${out_prefix}






