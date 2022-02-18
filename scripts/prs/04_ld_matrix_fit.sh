#!/usr/bin/env bash
#
#$ -N ld_matrix_fit
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/ld_matrix_fit.log
#$ -e logs/ld_matrix_fit.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 10
#$ -q short.qc@@short.hga
#$ -t 2
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/prs/04_ld_matrix_fit.R"

readonly bed_dir="data/prs/hapmap/ld"
readonly pheno_dir="data/phenotypes"
readonly sumstat_dir="data/prs/sumstat"
readonly out_dir="data/prs/hapmap/ld/corr"

readonly bed_file="${bed_dir}/short_merged_ukb_hapmap_rand_10k_eur.bed"
readonly pheno_file="${pheno_dir}/curated_phenotypes.tsv"

readonly pheno_list_cts="${pheno_dir}/curated_phenotypes_cts_header.tsv"
readonly phenotype_cts=$( cut -f${SGE_TASK_ID} ${pheno_list_cts} )
readonly pheno_list_binary="${pheno_dir}/curated_phenotypes_binary_header.tsv"
readonly phenotype_binary=$( cut -f${SGE_TASK_ID} ${pheno_list_binary} )

mkdir -p ${out_dir}

fit_ld_matrix() 
{
  SECONDS=0
  export OPENBLAS_NUM_THREADS=1 # avoid two levels of parallelization
  local sumstat_file="${sumstat_dir}/ukb_hapmap_500k_eur_${1}.txt.gz"
  local out_prefix="${out_dir}/ldpred2_${1}"
  set_up_rpy
  set -x
  Rscript "${rscript}" \
   --path_bed_ld "${bed_file}" \
   --path_sumstat ${sumstat_file} \
   --out_prefix ${out_prefix}
  set +x
  log_runtime ${SECONDS}
}

fit_ld_matrix ${phenotype_binary}



