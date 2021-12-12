#!/usr/bin/env bash
#
#$ -N merge_markers
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/submit_merge_markers.log
#$ -e logs/submit_merge_markers.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 37
#$ -V

source utils/vcf_utils.sh
source utils/hail_utils.sh

readonly mt_dir="derived/knockouts/211206"
readonly marker_dir="data/conditional/common/spa_conditional"
readonly out_dir="data/conditional/common/merge_markers"

readonly pheno_dir="data/phenotypes"
readonly pheno_list="${pheno_dir}/UKBB_WES200k_cts_phenotypes_header.txt"
readonly index=${SGE_TASK_ID}
readonly phenotype=$( cut -f${index} ${pheno_list} )

readonly merge_script="scripts/conditional/_merge_markers.sh"

mkdir -p ${out_dir}

submit_merge_markers_job() 
{
  if [ -f ${2} ]; then
    set -x
    qsub -N "_merge_markers}" \
      -t 9 \
      -q "short.qa" \
      -pe shmem 2 \
      "${merge_script}" \
      "${1}" \
      "${2}" \
      "${3}"
    set +x
  else
    echo >&2 "Conditioning markers file '${2}' does not exist."
  fi
}


annotation='ptv_damaging_missense'
mt="${mt_dir}/ukb_wes_200k_maf0_0.01_chrCHR_${annotation}_ko.mt"
markers_table="${marker_dir}/211111_spa_conditional_${annotation}_${phenotype}_conditioning_markers.txt"
out_prefix="${out_dir}/211111_${annotation}_${phenotype}_merged_chrCHR"
submit_merge_markers_job ${mt} ${markers_table} ${out_prefix}



