#!/usr/bin/env bash
#
#$ -N case_ctrl_icd_overview
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/case_ctrl_icd_overview.log
#$ -e logs/case_ctrl_icd_overview.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/14_case_ctrl_icd_overview.R"

readonly in_dir="data/post_hoc/results"
readonly in_file1="${in_dir}/ukb_eur_wes_200k_case_ctrl_pLoF_damaging_missense_by_phenotypes.txt.gz"
readonly in_file2="${in_dir}/ukb_eur_wes_200k_case_ctrl_pLoF_by_phenotypes.txt.gz"
readonly in_file3="${in_dir}/ukb_eur_wes_200k_case_ctrl_damaging_missense_by_phenotypes.txt.gz"

readonly icd_dir="data/phenotypes"
readonly icd="${icd_dir}/phenotype_icd_chapter.txt"

readonly out_dir="data/post_hoc/results"
readonly out_prefix="${out_dir}/case_ctrl_icd_overview"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --out_prefix "${out_prefix}" \
  --ko_by_phenotype1 "${in_file1}" \
  --ko_by_phenotype2 "${in_file2}" \
  --ko_by_phenotype3 "${in_file3}" \
  --icd_path "${icd}"

