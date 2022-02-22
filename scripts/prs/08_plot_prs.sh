#!/usr/bin/env bash
#
#$ -N plot_prs
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/plot_prs.log
#$ -e logs/plot_prs.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1
#$ -V

source utils/bash_utils.sh

readonly rscript="scripts/prs/08_plot_prs.R"

readonly prs_dir="data/prs/scores/test"
readonly prs="${prs_dir}/test_ukb_by_chrom_500k_pgs.txt.gz"

readonly pheno_dir="data/phenotypes"
readonly phenotypes="${pheno_dir}/curated_phenotypes_header.tsv"

readonly out_dir="data/prs/scores/test"
readonly out_prefix="${out_dir}/test_ukb_by_chrom_500k_pgs"

mkdir -p ${out_dir}

plot_prs()
{ 
  set_up_rpy
  set -x
  Rscript "${rscript}" \
   --phenotypes "${phenotypes}" \
   --prs "${prs}" \
   --out_prefix "${out_prefix}" 
  set +x
}

plot_prs




