#!/usr/bin/env bash
#
#$ -N ld_by_phenotypes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/ld_by_phenotypes.log
#$ -e logs/ld_by_phenotypes.errors.log
#$ -P lindgren.prjc
#$ -q short.qa
#$ -pe shmem 10
#$ -t 20-22
#$ -V

set -o errexit
set -o nounset



source utils/bash_utils.sh

calc_allele_count_by_phenotype() {
  set_up_rpy
  local out_prefix="${1}"
  local pheno_path="${2}"
  local pheno_file="${3}"
  local phenotype="${4}"
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





