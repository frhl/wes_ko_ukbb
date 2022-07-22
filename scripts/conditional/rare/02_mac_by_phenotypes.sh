#!/usr/bin/env bash
#
#$ -N mac_by_phenotypes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/mac_by_phenotypes.log
#$ -e logs/mac_by_phenotypes.errors.log
#$ -P lindgren.prjc
#$ -q short.qa
#$ -pe shmem 10
#$ -t 20-22
#$ -V

# Note: this script requires at least 10 A cores. For chromosome
# one, up to 20 A cores will be needed to run calc_allele_count_by_phenotype

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly rscript="scripts/conditional/rare/02_mac_by_phenotypes.R"
readonly rcombine="scripts/conditional/rare/_combine_ac.R"

readonly chr="${SGE_TASK_ID}"
readonly pheno_dir="data/phenotypes"
readonly in_dir="data/conditional/rare/combined"
readonly out_dir="data/conditional/rare/combined"

readonly pheno_cts_path="${pheno_dir}/curated_covar_phenotypes_cts_200k.tsv"
readonly pheno_bin_path="${pheno_dir}/curated_covar_phenotypes_binary_200k.tsv"

readonly input_path="${in_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense.vcf.bgz"
readonly out_prefix_cts="${out_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense_AC_cts"
readonly out_prefix_bin="${out_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense_AC_bin"
readonly out_prefix_combined="${out_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense_AC"
readonly covar_path="${pheno_dir}/covars1.csv"

readonly phenotypes_cts=$(cat "${pheno_dir}/filtered_phenotypes_cts_manual.tsv" | tr "\n" "," | sed 's/\(.*\),/\1 /' )
readonly phenotypes_bin=$(cat "${pheno_dir}/filtered_phenotypes_binary_header.tsv" | tr "\n" "," | sed 's/\(.*\),/\1 /' )

set_up_rpy

calc_allele_count_by_phenotype() {
  local out_prefix="${1}"
  local pheno_path="${2}"
  local phenotypes="${3}"
  if [ ! -f "${out_prefix}.txt.gz" ]; then
  Rscript "${rscript}" \
     --in_vcf ${input_path} \
     --out_prefix ${out_prefix} \
     --subset_phenotypes ${phenotypes} \
     --phenotypes ${pheno_path} \
     --covariates ${covar_path}
  else
    >&2 echo "${out_prefix}.txt.gz already exists. Skipping."
  fi
}


# calc AC for binary/cts 
calc_allele_count_by_phenotype ${out_prefix_cts} ${pheno_cts_path} ${phenotypes_cts}
calc_allele_count_by_phenotype ${out_prefix_bin} ${pheno_bin_path} ${phenotypes_bin}

# combine the two files
Rscript "${rcombine}" \
  --file_bin "${out_prefix_bin}.txt.gz" \
  --file_cts "${out_prefix_cts}.txt.gz" \
  --out_prefix "${out_prefix_combined}"





