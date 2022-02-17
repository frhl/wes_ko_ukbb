#!/usr/bin/env bash
#
#$ -N prs
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prs.log
#$ -e logs/prs.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 2
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly bash_script="scripts/prs/_prs.sh"
readonly rscript="scripts/prs/05_prs.R"

readonly pred_dir="data/prs/hapmap"
readonly ld_dir="data/prs/hapmap/ld"
readonly sumstat_dir="data/prs/sumstat"
readonly pheno_dir="data/phenotypes"
readonly out_dir="data/prs/scores"

readonly pheno_list_cts="${pheno_dir}/curated_phenotypes_cts_header.tsv"
readonly phenotype_cts=$( cut -f${SGE_TASK_ID} ${pheno_list_cts} )
readonly pheno_list_binary="${pheno_dir}/curated_phenotypes_binary_header.tsv"
readonly phenotype_binary=$( cut -f${SGE_TASK_ID} ${pheno_list_binary} )

mkdir -p ${out_dir}

calc_prs_by_chrom()
{ 
  export OPENBLAS_NUM_THREADS=1 # avoid two levels of parallelization
  local phenotype=${1}
  local pred="${pred_dir}/ukb_hapmap_500k_eur_chrCHR.bed"
  local ld="${ld_dir}/short_ukb_hapmap_rand_10k_eur_chrCHR.bed"
  local sumstat="${sumstat_dir}/ukb_hapmap_500k_eur_${phenotype}_chrCHR.txt.gz"
  local out_prefix="${out_dir}/ukb_eur_prs_${phenotype}_chrCHR"
  set -x
  qsub -N "_prs_${phenotype}" \
    -t 21 \
    -q short.qc@@short.hge \
    -pe shmem 1 \
    "${bash_script}" \
    "${rscript}" \
    "${ld}" \
    "${sumstat}" \
    "${pred}" \
    "${out_prefix}"
  set +x 
}

calc_prs_by_chrom ${phenotype_binary}

