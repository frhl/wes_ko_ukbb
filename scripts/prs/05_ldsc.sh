#!/usr/bin/env bash
#
#$ -N ldsc
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/ldsc.log
#$ -e logs/ldsc.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/prs/05_ldsc.R"

readonly gwas_dir="data/prs/sumstat"
readonly bed_dir="data/prs/hapmap/ld/unrel_eur_10k"
readonly out_dir="data/prs/ldsc"
readonly pheno_dir="data/phenotypes"

readonly ld_bed="${bed_dir}/short_merged_ukb_hapmap_rand_10k_eur.bed"
readonly ld_dir="data/prs/hapmap/ld/matrix"

readonly index=${SGE_TASK_ID}

readonly file_cts="${pheno_dir}/filtered_phenotypes_cts.tsv"
readonly pheno_list_cts="${pheno_dir}/filtered_phenotypes_cts_header.tsv"
readonly phenotype_cts=$( sed "${index}q;d" ${pheno_list_cts} )

readonly file_binary="${pheno_dir}/filtered_phenotypes_binary.tsv"
readonly pheno_list_binary="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
readonly phenotype_binary=$( sed "${index}q;d" ${pheno_list_binary} )


mkdir -p ${out_dir}

export OPENBLAS_NUM_THREADS=1 # avoid two levels of parallelization

estimate_heritability(){
  set_up_rpy
  if [ ! -f "${out_prefix}.rds" ]; then
    local phenotype="${1}" 
    local trait="${2}"
    local out_prefix="${out_dir}/ldsc_${phenotype}"
    local gwas="${gwas_dir}/${trait}/ukb_hapmap_500k_eur_${phenotype}.txt.gz"
    set -x
    Rscript "${rscript}" \
        --gwas "${gwas}" \
        --ld_bed "${ld_bed}" \
        --ld_dir "${ld_dir}" \
        --trait "${trait}" \
        --out_prefix "${out_prefix}"
    set +x
  else
    echo "Note: ${out_prefix} already exists. Skipping.."
  fi
}

estimate_heritability "${phenotype_binary}" "binary"
#estimate_heritability "${phenotype_cts}" "cts"

