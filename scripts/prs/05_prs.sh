#!/usr/bin/env bash
#
#$ -N prs
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prs.log
#$ -e logs/prs.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc@@short.hge
#$ -t 1
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/prs/04_ld_matrix_fit.R"

readonly in_dir="data/prs/hapmap/ld"
readonly ld_dir="data/prs/hapmap/ld"
readonly sumstat_dir="data/prs/sumstat"
readonly pheno_dir="data/phenotypes"
readonly out_dir="data/prs/hapmad/ld"

readonly ld_matrix="${ld_dir}/"
readonly pheno_file="${pheno_dir}/curated_phenotypes.tsv"

readonly pheno_list_cts="${pheno_dir}/curated_phenotypes_cts_header.tsv"
readonly phenotype_cts=$( cut -f${SGE_TASK_ID} ${pheno_list_cts} )
readonly pheno_list_binary="${pheno_dir}/curated_phenotypes_binary_header.tsv"
readonly phenotype_binary=$( cut -f${SGE_TASK_ID} ${pheno_list_binary} )

readonly out_prefix="${out_dir}/ukb_eur_ld_10k_${phenotype}"

mkdir -p ${out_dir}

get_prs_matrix()
{
  local sumstat="${sumstat_dir}/ukb_hapmap_500k_eur_${1}_combined.txt.gz"
  set_up_rpy
  set -x
  Rscript "${rscript}" \
   --path_bed_ld "${bed}" \
   --path_sumstat ${sumstat} \
   --out_prefix ${out_prefix}
  set +x
}

fit_ld_matrix ${phenotype_binary}

