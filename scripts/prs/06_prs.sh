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
readonly rscript="scripts/prs/06_prs.R"

readonly ldsc_dir="data/prs/ldsc"
readonly pred_dir="data/prs/hapmap/ukb_500k"
readonly ld_dir="data/prs/hapmap/ld/matrix"
readonly pheno_dir="data/phenotypes"
readonly out_dir="data/prs/scores/test"

readonly pheno_list_cts="${pheno_dir}/curated_phenotypes_cts_header.tsv"
readonly phenotype_cts=$( cut -f${SGE_TASK_ID} ${pheno_list_cts} )
readonly pheno_list_binary="${pheno_dir}/curated_phenotypes_binary_header.tsv"
readonly phenotype_binary=$( cut -f${SGE_TASK_ID} ${pheno_list_binary} )

readonly method="inf"

mkdir -p ${out_dir}

export OPENBLAS_NUM_THREADS=1 # avoid two levels of parallelization

calc_prs_by_chrom()
{ 
  local phenotype=${1}
  local pred="${pred_dir}/ukb_hapmap_500k_eur_chrCHR.bed"
  local ldsc="${ldsc_dir}/ldsc_${phenotype}.rds"
  local out_prefix="${out_dir}/prs_inf_${phenotype}_chrCHR"
  set -x
  qsub -N "_prs_${phenotype}" \
    -t 21 \
    -q short.qc@@short.hga \
    -pe shmem 1 \
    "${bash_script}" \
    "${rscript}" \
    "${pred}" \
    "${ldsc}" \
    "${ld_dir}" \
    "${method}" \
    "${out_prefix}"
  set +x 
}

calc_prs_by_chrom ${phenotype_binary}

