#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=ldsc
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/ldsc.log
#SBATCH --error=logs/ldsc.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 4
#SBATCH --array=1-50

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/prs/05_ldsc.R"

readonly gwas_dir="data/prs/sumstat"
readonly bed_dir="data/prs/hapmap/ld/unrel_kin_eur_10k"
readonly out_dir="data/prs/ldsc_test"
readonly pheno_dir="data/phenotypes"

readonly ld_bed="${bed_dir}/short_merged_ukb_hapmap_rand_10k_eur.bed"
readonly ld_dir="data/prs/hapmap/ld/matrix_unrel_kin"

readonly cluster=$( get_current_cluster )
readonly index=$( get_array_task_id )

readonly file_cts="${pheno_dir}/curated_covar_phenotypes_cts_int_500k.txt.gz"
readonly pheno_list_cts="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
readonly phenotype_cts=$( sed "${index}q;d" ${pheno_list_cts} )

readonly file_binary="${pheno_dir}/dec22_phenotypes_binary_500k.tsv.gz"
readonly pheno_list_binary="${pheno_dir}/dec22_phenotypes_binary_200k_header.tsv"
readonly phenotype_binary=$( sed "${index}q;d" ${pheno_list_binary} )

mkdir -p ${out_dir}

export OPENBLAS_NUM_THREADS=1 # avoid two levels of parallelization

estimate_heritability(){
  set_up_ldpred2
  local phenotype="${1}" 
  local trait="${2}"
  local out_prefix="${out_dir}/ldsc_${phenotype}"
  if [ ! -f "${out_prefix}.rds" ]; then
    local gwas="${gwas_dir}/${trait}/ukb_hapmap_500k_eur_${phenotype}.txt.gz"
    set -x
    Rscript "${rscript}" \
        --gwas "${gwas}" \
        --ld_bed "${ld_bed}" \
        --ld_dir "${ld_dir}" \
        --trait "${trait}" \
        --phenotype "${phenotype}" \
        --path_cts_phenotypes "${file_cts}" \
        --out_prefix "${out_prefix}" \
        --disable_qc # we have already QCed Sumstats/SNPs
    set +x
  else
    echo "Note: ${out_prefix} already exists. Skipping.."
  fi
}


estimate_heritability "${phenotype_cts}" "cts"
#estimate_heritability "${phenotype_binary}" "binary"

