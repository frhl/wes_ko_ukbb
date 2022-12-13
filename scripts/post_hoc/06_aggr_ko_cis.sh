#!/usr/bin/env bash
#
#$ -N aggr_ko_cis
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/aggr_ko_cis.log
#$ -e logs/aggr_ko_cis.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/06_aggr_ko_cis.R"

readonly knockout_dir="data/knockouts/alt/pp90/combined"

readonly out_dir="data/knockouts/tables"
readonly out_prefix="${out_dir}/221312_aggr_ko_table"

mkdir -p ${out_dir}


aggr_knockouts() {
  local annotation=${1}
  local pattern=${2}
  local out_prefix_annotation="${out_prefix}_${annotation}"
  set_up_rpy
  Rscript "${rscript}" \
     --knockout_dir ${knockout_dir} \
     --knockout_pattern ${pattern} \
     --out_prefix "${out_prefix_annotation}"
}

aggr_knockouts "pLoF_damaging_missense" "chr[0-9]+_pLoF_damaging_missense_all.tsv.gz"
aggr_knockouts "pLoF" "chr[0-9]+_pLoF_all.tsv.gz"
aggr_knockouts "damaging_missense" "chr[0-9]+_damaging_missense_all.tsv.gz"

