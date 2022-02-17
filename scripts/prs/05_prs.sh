#!/usr/bin/env bash
#
#$ -N prs
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prs.log
#$ -e logs/prs.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 1
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly bash_script="scripts/prs/_prs.sh"
readonly rscript="scripts/prs/05_prs.R"

readonly pred_dir="data/prs/hapmap"
readonly ld_matrix_dir="data/prs/hapmap/ld"
readonly sumstat_dir="data/prs/sumstat"
readonly pheno_dir="data/phenotypes"
readonly out_dir="data/prs/hapmad/ld"

readonly pred="${pred_dir}/ukb_hapmap_500k_eur_chrCHR.bed"

readonly pheno_list_cts="${pheno_dir}/curated_phenotypes_cts_header.tsv"
readonly phenotype_cts=$( cut -f${SGE_TASK_ID} ${pheno_list_cts} )
readonly pheno_list_binary="${pheno_dir}/curated_phenotypes_binary_header.tsv"
readonly phenotype_binary=$( cut -f${SGE_TASK_ID} ${pheno_list_binary} )

mkdir -p ${out_dir}

calc_prs_by_chrom()
{ 
  local ld_matrix="${ld_matrix_dir}/ukb_eur_ld_10k_${1}.rds"
  local sumstat="${sumstat_dir}/ukb_hapmap_500k_eur_${1}_combined.txt.gz"
  local out_prefix="${out_dir}/ukb_eur_prs_${1}_chrCHR"
  set -x
  qsub -N "_prs_${phenotype}" \
    -t 21 \
    -q short.qc@@short.hge \
    -pe shmem 1 \
    "${bash_script}" \
    "${rscript}" \
    "${ld_matrix}" \
    "${sumstat}" \
    "${pred}" \
    "${out_prefix}"
  set +x 
}

calc_prs_by_chrom ${phenotype_binary}

