#!/usr/bin/env bash
#
#$ -N merge_markers
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/merge_markers.log
#$ -e logs/merge_markers.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q short.qc
#$ -t 3

source utils/vcf_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly marker_dir=""
readonly vcf_dir=""
readonly out_dir="data/conditional/common"

readonly hail_script="scripts/conditional/02_filter_genotypes.py"
readonly vcf="${vcf_dir}/"
readonly markers="${marker_dir}/"

readonly out_prefix="${out_dir}/"

set_up_hail
set_up_pythonpath_legacy
mkdir -p ${out_dir}

merge_markers() {
  SECONDS=0
  local vcf=${1}
  local markers=${2}
  local out_prefix=${3}
  python3 "${hail_script}" \
     --vcf ${vcf} \
     --markers ${markers} \
     --out_prefix ${out_prefix} \ 
  print_update "Hail finished writing VCFs in ${SECONDS}"
}


annotation='synonymous'
out_prefix="${out_dir}/211111_intervals_${annotation}_${phenotype}"
merge_markers ${annotation} ${out_prefix}



