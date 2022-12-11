#!/usr/bin/env bash
#
#$ -N logit_burden
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/logit_burden.log
#$ -e logs/logit_burden.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 3
#$ -q short.qc
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/05_logit_burden.R"

readonly unrelated_dir="data/post_hoc/unrelated"
readonly unrelated_path="${unrelated_dir}/ukb_wes_ko_samples.txt.gz"

readonly phenotypes_dir="data/phenotypes"
readonly phenotypes_path="${phenotypes_dir}/spiros_brava_phenotypes_binary_200k.tsv.gz"
readonly phenotypes_header="${phenotypes_dir}/spiros_brava_phenotypes_binary_200k_header.tsv"
readonly covars_path="${phenotypes_dir}/covars1.csv"

readonly knockout_dir="data/knockouts/alt/pp90/combined"
readonly knockout_regex="pLoF_damaging_missense_all"

readonly out_dir="data/knockouts/tables"
readonly out_prefix="${out_dir}/pLoF_damaging_missense_full_burden_logit"

mkdir -p ${out_dir}


calc_burden() {
  local variable=${1}
  local out_prefix_variable="${out_prefix}.${variable}"
  set_up_rpy
  Rscript "${rscript}" \
   --phenotypes_path ${phenotypes_path} \
   --phenotypes_header ${phenotypes_header} \
   --unrelated_path ${unrelated_path} \
   --covars_path ${covars_path} \
   --knockout_dir ${knockout_dir} \
   --knockout_pattern ${knockout_regex} \
   --out_prefix "${out_prefix_variable}" \
   --variable "${variable}"
}

calc_burden "het"
calc_burden "ko"
calc_burden "chet"
calc_burden "homs"


