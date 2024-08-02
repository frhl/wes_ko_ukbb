#!/usr/bin/env bash
#
#$ -N ldsc_auc_icd
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/ldsc_auc_icd.log
#$ -e logs/ldsc_auc_icd.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/06_ldsc_auc_icd.R"


readonly summary_ldsc="data/prs/validation/ldsc_summary.txt.gz"
readonly summary_prs_bin="data/prs/validation/spiro_pgs_auc_summary.txt.gz"

readonly icd_dir="data/phenotypes"
readonly icd="${icd_dir}/phenotype_icd_chapter.txt"

readonly out_dir="data/post_hoc/results"
readonly out_prefix="${out_dir}/ldsc_auc_icd"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --out_prefix "${out_prefix}" \
  --summary_ldsc "${summary_ldsc}" \
  --summary_prs_bin "${summary_prs_bin}" \
  --icd "${icd}"

