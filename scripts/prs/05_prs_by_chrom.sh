#!/usr/bin/env bash
#
#$ -N prs_by_chrom
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prs_by_chrom.log
#$ -e logs/prs_by_chrom.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 2
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly bash_script="scripts/prs/_prs_by_chrom.sh"
readonly rscript="scripts/prs/05_prs_by_chrom.R"

readonly pred_dir="data/prs/hapmap/ukb_500k"
readonly gwas_dir="data/prs/sumstat"
readonly bed_dir="data/prs/hapmap/ld/unrel_eur_10k"
readonly ld_dir="data/prs/hapmap/ld/matrix"
readonly pheno_dir="data/phenotypes"
readonly out_dir="data/prs/scores/by_chrom"

readonly ld_bed="${bed_dir}/short_merged_ukb_hapmap_rand_10k_eur.bed"

# setup paths to phenotypes
readonly pheno_list_cts="${pheno_dir}/curated_phenotypes_cts_header.tsv"
readonly phenotype_cts=$( cut -f${SGE_TASK_ID} ${pheno_list_cts} )
readonly pheno_list_binary="${pheno_dir}/curated_phenotypes_binary_header.tsv"
readonly phenotype_binary=$( cut -f${SGE_TASK_ID} ${pheno_list_binary} )

# what ldpred2 method should be used ("auto" or "inf")
readonly method="inf"

# what trait is being considered? This affects how n_eff is 
# calculated. Should be either "binary" or "cts"
readonly trait="binary"

mkdir -p ${out_dir}

export OPENBLAS_NUM_THREADS=1 # avoid two levels of parallelization

calc_prs_by_chrom()
{ 
  local phenotype=${1}
  local pred="${pred_dir}/ukb_hapmap_500k_eur_chrCHR.bed"
  local gwas="${gwas_dir}/ukb_hapmap_500k_eur_${phenotype}.txt.gz"
  local out_prefix="${out_dir}/prs_inf_${phenotype}_chrCHR"
  set -x
  qsub -N "_prs_c_${phenotype}" \
    -t 21 \
    -q short.qc@@short.hga \
    -pe shmem 2 \
    "${bash_script}" \
    "${rscript}" \
    "${gwas}" \
    "${pred}" \
    "${ld_bed}" \
    "${ld_dir}" \
    "${method}" \
    "${trait}" \
    "${out_prefix}"
  set +x 
}

calc_prs_by_chrom ${phenotype_binary}

