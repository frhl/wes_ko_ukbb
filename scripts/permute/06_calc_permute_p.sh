#!/usr/bin/env bash
#
#$ -N calc_permute_p
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/calc_permute_p.log
#$ -e logs/calc_permute_p.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 44
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"
readonly bash_script="scripts/permute/_array_spa.sh"

readonly alt_dir="data/saige/output/....."
readonly null_dir="data/permute/permutations/spa/CHR"
readonly out_dir="data/permute/empirical_p/"


readonly in_prefix="${in_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_permuted_chrCHR_GENE"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_permuted_chrCHR_GENE"

readonly overview="${overview_dir}/overview.tsv.gz"
readonly gene_pheno_p="${overview_dir}/overview_pvalue_genes.tsv.gz"
readonly gene_spa="${overview_dir}/overview_genes.tsv.gz"
readonly annotation="pLoF_damaging_missense"

calc_p_binary()
{
  local pheno_list="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  calc_p "${phenotype}" "binary"
}

calc_p_cts()
{
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  calc_p  "${phenotype}" "cts"
}

calc_p()
{

  local phenotype=${1?Error: Missing arg2 (phenotype)}
  local trait=${2?Error: Missing arg3 (trait)}
  local out_prefix_pheno="${out_prefix}_${phenotype}"
  mkdir -p ${step2_dir}
  set -x
  qsub -N "spa_${phenotype}" \
    -t ${tasks} \
    -q "test.qc" \
    -pe shmem ${nslots} \
    "${bash_script}" \
    "${phenotype}" \
    "${in_prefix}" \
    "${in_gmat}" \
    "${in_var}" \
    "${min_mac}" \
    "${annotation}" \
    "${gene_pheno_p}" \
    "${overview}" \
    "${gene_spa}" \
    "${p_per_job}" \
    "${out_prefix_pheno}" \
    "${out_dir}"
}


