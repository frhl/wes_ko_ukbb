#!/usr/bin/env bash


# path to bonferonni corrected phenotypes
get_prs_path() {
  echo "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/prs/validation/ldsc_summary_bonf_sig_phenos.txt"
}

# get path to currently used phenotype headers and phenotypes
get_pheno_header_path() {
 echo "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/phenotype/dec22_phenotypes_binary_200k_header.tsv"
}

get_pheno_header_path_200k() {
  echo "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/phenotypes/dec22_phenotypes_binary_200k.tsv.gz"
}

# check if we allow conditioning om prs
pheno_allow_cond_prs() {
  local phenotype="${1}"
  echo "$( cat $(get_prs_path) | grep -w "${phenotype}" | wc -l)"
}

# check if phenotype 
validate_phenotype() {
  local phenotype="${1}"
  local ok="$( cat $(get_pheno_header_path) | grep -w "${phenotype}" | wc -l)"
  if [ "${ok}" -eq "0" ]; then
    raise_error "${phenotype} not pheno path!"
  fi
}











